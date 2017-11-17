let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "ac355040656de04f59406ba2380a96f4124ebdad";
    sha256 = "0frhc7mnx88sird6ipp6578k5badibsl0jfa22ab9w6qrb88j825";
  });
in
  with import nixpkgs {
    overlays = [(self: super:
    {
      atlas = super.atlas.override {
        withLapack = true;
      };
      hdf5 = super.hdf5.override {
        gfortran = super.gfortran;
        mpi = super.openmpi;
      };
    }
    )];
  };
  stdenv.mkDerivation {
    name = "HANDE";
    buildInputs = [
      atlas
      clang-tools
      cmake
      exa
      gcc
      gdb
      gfortran
      hdf5
      libuuid
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
      valgrind
    ];
  }
