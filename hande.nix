with import <nixpkgs> {}; {
  handeEnv = stdenv.mkDerivation {
    name = "HANDE";
    buildInputs = [
      atlasWithLapack
      ccache
      clang
      clang-analyzer
      clang-tools
      cmake
      gcc
      gdb
      gfortran
      hdf5
      hdf5-mpi
      libuuid
      lldb
      lua5_3
      openmpi
      python35Packages.jupyter
      python35Packages.matplotlib
      python35Packages.numpy
      python35Packages.pandas
      python35Packages.pyyaml
      python35Packages.recommonmark
      python35Packages.scipy
      python35Packages.sphinx
      python35Packages.sphinx_rtd_theme
      python35Packages.sympy
    ];
  };
}
