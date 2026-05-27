#include "src/Bin.h"
#include "src/BwhamBinning.h"
#include "src/SimpleBias.h"
#include "src/SquaredBias.h"
#include "tools/CommandLineArguments.h"
#include "tools/InputParser.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
bool near(double actual, double expected, double tolerance)
{
    return std::abs(actual - expected) <= tolerance;
}

int require_true(const std::string& name, bool condition)
{
    if (!condition)
    {
        std::cerr << name << " expected true\n";
        return 1;
    }

    return 0;
}

int require_false(const std::string& name, bool condition)
{
    if (condition)
    {
        std::cerr << name << " expected false\n";
        return 1;
    }

    return 0;
}

template <typename T>
int require_equal(const std::string& name, const T& actual, const T& expected)
{
    if (actual != expected)
    {
        std::cerr << name << " expected " << expected << " but got " << actual << "\n";
        return 1;
    }

    return 0;
}

int require_near(const std::string& name, double actual, double expected, double tolerance)
{
    if (!near(actual, expected, tolerance))
    {
        std::cerr << name << " expected " << expected << " but got " << actual << "\n";
        return 1;
    }

    return 0;
}

ParameterPack make_bin_pack()
{
    ParameterPack pack("bins");
    pack.insert("range", std::vector<std::string>{"0.0", "10.0"});
    pack.insert("numbins", "5");
    pack.insert("dimension", "2");
    return pack;
}

ParameterPack make_simple_bias_pack()
{
    ParameterPack pack("bias");
    pack.insert("dimension", "2");
    pack.insert("kappa", std::vector<std::string>{"2.0", "4.0"});
    pack.insert("xstar", std::vector<std::string>{"1.0", "-1.0"});
    pack.insert("phi", std::vector<std::string>{"0.5", "-0.25"});
    return pack;
}

ParameterPack make_squared_bias_pack()
{
    ParameterPack pack("bias");
    pack.insert("dimension", "2");
    pack.insert("phi", std::vector<std::string>{"2.0", "0.5"});
    return pack;
}

ParameterPack make_bin_pack_with_dimension(
    const std::vector<std::string>& range,
    const std::string& numbins,
    const std::string& dimension)
{
    ParameterPack pack("bins");
    pack.insert("range", range);
    pack.insert("numbins", numbins);
    pack.insert("dimension", dimension);
    return pack;
}

int test_input_parser()
{
    const std::string input_name = "test_core_utilities_input.dat";
    {
        std::ofstream input(input_name);
        input << "# comment line\n";
        input << "title = parser_smoke\n";
        input << "temperature = 310.5\n";
        input << "enabled = true\n";
        input << "values = [ 1.0 -2.5 3.25 ]\n";
        input << "timeseries = {\n";
        input << "  path = first.dat\n";
        input << "  columns = [ 1 3 ]\n";
        input << "}\n";
        input << "timeseries = {\n";
        input << "  path = second.dat\n";
        input << "  columns = [ 2 ]\n";
        input << "}\n";
    }

    ParameterPack pack("root");
    InputParser parser;
    parser.ParseFile(input_name, pack);
    std::remove(input_name.c_str());

    int failures = 0;
    std::string title;
    double temperature = 0.0;
    bool enabled = false;
    std::vector<double> values;

    failures += require_true(
        "InputParser reads string",
        pack.ReadString("title", ParameterPack::KeyType::Required, title));
    failures += require_equal("InputParser string value", title, std::string("parser_smoke"));

    failures += require_true(
        "InputParser reads number",
        pack.ReadNumber("temperature", ParameterPack::KeyType::Required, temperature));
    failures += require_near("InputParser number value", temperature, 310.5, 1e-12);

    failures += require_true(
        "InputParser reads bool",
        pack.Readbool("enabled", ParameterPack::KeyType::Required, enabled));
    failures += require_true("InputParser bool value", enabled);

    failures += require_true(
        "InputParser reads numeric vector",
        pack.ReadVectorNumber("values", ParameterPack::KeyType::Required, values));
    failures += require_equal("InputParser vector size", static_cast<int>(values.size()), 3);
    failures += require_near("InputParser vector first", values[0], 1.0, 1e-12);
    failures += require_near("InputParser vector second", values[1], -2.5, 1e-12);
    failures += require_near("InputParser vector third", values[2], 3.25, 1e-12);

    std::string missing;
    failures += require_false(
        "InputParser optional missing string",
        pack.ReadString("missing", ParameterPack::KeyType::Optional, missing));

    auto timeseries = pack.findParamPacks("timeseries", ParameterPack::KeyType::Required);
    failures += require_equal("InputParser repeated pack count", static_cast<int>(timeseries.size()), 2);

    std::string first_path;
    std::vector<int> second_columns;
    failures += require_true(
        "InputParser first nested path",
        timeseries[0]->ReadString("path", ParameterPack::KeyType::Required, first_path));
    failures += require_equal("InputParser first nested path value", first_path, std::string("first.dat"));
    failures += require_true(
        "InputParser second nested columns",
        timeseries[1]->ReadVectorNumber("columns", ParameterPack::KeyType::Required, second_columns));
    failures += require_equal("InputParser second nested column count", static_cast<int>(second_columns.size()), 1);
    failures += require_equal("InputParser second nested column value", second_columns[0], 2);

    return failures;
}

int test_command_line_arguments()
{
    const char* raw_argv[] = {
        "Wham",
        "input.dat",
        "-abspath",
        "testdata",
        "-window",
        "-1.5",
        "-0.25",
        "--label",
        "sample"
    };
    std::vector<char*> argv;
    for (const char* arg : raw_argv)
    {
        argv.push_back(const_cast<char*>(arg));
    }

    CommandLineArguments args(static_cast<int>(argv.size()), argv.data());

    int failures = 0;
    std::string abspath;
    std::string label;
    std::vector<double> window;

    failures += require_equal("CommandLineArguments key count", args.get_num_keys(), 3);
    failures += require_true("CommandLineArguments has abspath", args.has_key("abspath"));
    failures += require_true("CommandLineArguments has label", args.has_key("label"));
    failures += require_false("CommandLineArguments missing key", args.has_key("missing"));
    failures += require_true(
        "CommandLineArguments reads string",
        args.readString("abspath", CommandLineArguments::Keys::Required, abspath));
    failures += require_equal("CommandLineArguments string value", abspath, std::string("testdata"));
    failures += require_true(
        "CommandLineArguments reads double vector",
        args.readVector("window", CommandLineArguments::Keys::Required, window));
    failures += require_equal("CommandLineArguments vector size", static_cast<int>(window.size()), 2);
    failures += require_near("CommandLineArguments negative first", window[0], -1.5, 1e-12);
    failures += require_near("CommandLineArguments negative second", window[1], -0.25, 1e-12);
    failures += require_true(
        "CommandLineArguments double dash key",
        args.readString("label", CommandLineArguments::Keys::Required, label));
    failures += require_equal("CommandLineArguments double dash value", label, std::string("sample"));

    return failures;
}

int test_bin()
{
    Bin bin(make_bin_pack());

    int failures = 0;
    failures += require_equal("Bin dimension", bin.getDimension(), 2);
    failures += require_equal("Bin count", bin.getNumbins(), 5);
    failures += require_near("Bin step", bin.getStep(), 2.0, 1e-12);
    failures += require_true("Bin includes lower bound", bin.isInRange(0.0));
    failures += require_true("Bin includes interior", bin.isInRange(9.999));
    failures += require_false("Bin excludes upper bound", bin.isInRange(10.0));
    failures += require_false("Bin excludes below range", bin.isInRange(-0.001));
    failures += require_equal("Bin lower index", bin.findBin(0.0), 0);
    failures += require_equal("Bin interior index", bin.findBin(2.0), 1);
    failures += require_equal("Bin last index", bin.findBin(9.999), 4);
    failures += require_near("Bin center", bin.getLocationOfBin(2), 5.0, 1e-12);

    return failures;
}

int test_bwham_binning_helpers()
{
    ParameterPack dim1_pack = make_bin_pack_with_dimension({"0.0", "2.0"}, "2", "1");
    ParameterPack dim2_pack = make_bin_pack_with_dimension({"-1.0", "1.0"}, "4", "2");
    Bin dim1(dim1_pack);
    Bin dim2(dim2_pack);
    std::vector<const Bin*> bins{&dim1, &dim2};

    BwhamBinning::BinGrid grid = BwhamBinning::buildBinGrid(bins, 2);

    int failures = 0;
    failures += require_equal("BwhamBinning total bins", grid.totalBins, 8);
    failures += require_equal("BwhamBinning center count", static_cast<int>(grid.centers.size()), 8);

    std::vector<int> lower_index{0, 0};
    std::vector<int> upper_index{1, 3};
    failures += require_equal("BwhamBinning lower flat index", grid.indexToFlat.at(lower_index), 0);
    failures += require_equal("BwhamBinning upper flat index", grid.indexToFlat.at(upper_index), 7);
    failures += require_near("BwhamBinning lower center dim1", grid.centers[0][0], 0.5, 1e-12);
    failures += require_near("BwhamBinning lower center dim2", grid.centers[0][1], -0.75, 1e-12);
    failures += require_near("BwhamBinning upper center dim1", grid.centers[7][0], 1.5, 1e-12);
    failures += require_near("BwhamBinning upper center dim2", grid.centers[7][1], 0.75, 1e-12);

    std::vector<int> sample_index;
    failures += require_true(
        "BwhamBinning maps sample by declared dimensions",
        BwhamBinning::findBinIndexForSample(bins, {1.2, 0.6}, sample_index));
    failures += require_equal("BwhamBinning sample dim1 index", sample_index[0], 1);
    failures += require_equal("BwhamBinning sample dim2 index", sample_index[1], 3);

    std::vector<const Bin*> reversed_bins{&dim2, &dim1};
    failures += require_true(
        "BwhamBinning maps sample independent of bin order",
        BwhamBinning::findBinIndexForSample(reversed_bins, {1.2, 0.6}, sample_index));
    failures += require_equal("BwhamBinning reversed sample dim1 index", sample_index[0], 1);
    failures += require_equal("BwhamBinning reversed sample dim2 index", sample_index[1], 3);

    failures += require_false(
        "BwhamBinning rejects upper bound",
        BwhamBinning::findBinIndexForSample(bins, {1.2, 1.0}, sample_index));

    return failures;
}

int test_biases()
{
    SimpleBias simple(make_simple_bias_pack());
    SquaredBias squared(make_squared_bias_pack());

    int failures = 0;

    const std::vector<double> simple_x{3.0, 2.0, 9.0};
    std::vector<double> simple_force = simple.calculateForce(simple_x);
    failures += require_near("SimpleBias energy", simple.calculate(simple_x), 23.0, 1e-12);
    failures += require_equal("SimpleBias force size", static_cast<int>(simple_force.size()), 3);
    failures += require_near("SimpleBias force first", simple_force[0], -4.5, 1e-12);
    failures += require_near("SimpleBias force second", simple_force[1], -11.75, 1e-12);
    failures += require_near("SimpleBias force untouched dimension", simple_force[2], 0.0, 1e-12);

    const std::vector<double> squared_x{3.0, -4.0, 5.0};
    std::vector<double> squared_force = squared.calculateForce(squared_x);
    failures += require_near("SquaredBias energy", squared.calculate(squared_x), 26.0, 1e-12);
    failures += require_equal("SquaredBias force size", static_cast<int>(squared_force.size()), 3);
    failures += require_near("SquaredBias force first", squared_force[0], -12.0, 1e-12);
    failures += require_near("SquaredBias force second", squared_force[1], 4.0, 1e-12);
    failures += require_near("SquaredBias force untouched dimension", squared_force[2], 0.0, 1e-12);

    return failures;
}
}

int main()
{
    int failures = 0;
    failures += test_input_parser();
    failures += test_command_line_arguments();
    failures += test_bin();
    failures += test_bwham_binning_helpers();
    failures += test_biases();

    return failures == 0 ? 0 : 1;
}
