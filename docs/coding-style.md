# Coding Style

This project uses `clang-format` for C++ formatting. Run it from the repository root
before committing source changes:

```bash
find src tools parallel test -name '*.cpp' -o -name '*.h' | sort | xargs clang-format -i
```

Naming conventions:

- Classes, structs, and type aliases use `PascalCase`.
- Functions and methods use `camelCase`.
- Local variables and function parameters use `camelCase`.
- Private and protected data members use `camelCase_`.
- Input-file keywords, registered output names, and scientific symbols may keep their
  established spelling when renaming would break user-facing behavior or obscure the
  equation being implemented.

Avoid formatting vendored dependencies such as `Eigen/` and `LBFGS/`.
