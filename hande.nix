with import <nixpkgs> {}; {
  handeEnv = stdenv.mkDerivation {
    name = "HANDE";
    buildInputs = [
      ccache
      clang
      clang-analyzer
      clang-tools
      cmake
      gcc
      gdb
      gfortran
      hdf5
      libuuid
      lldb
      lua5_3
      openmpi
      python35Packages.jupyter
      python35Packages.matplotlib
      python35Packages.numpy
      python35Packages.pandas
      python35Packages.scipy
      python35Packages.sphinx
      python35Packages.sphinx_rtd_theme
      python35Packages.sympy
    ];
  };
}
