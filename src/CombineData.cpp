#include "CombineData.h"

namespace timeseriesOP
{
    registry<CombineData> registerCD("combine_data");
}

CombineData::CombineData(const TSInput& input)
: TSoperation(input)
{

}