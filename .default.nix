let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-unstable";
    sha256 = "0r92p2blm9k2dhy5xv41k9asqf9xfd6if4c801ab1wjzy8wsf1wi";
  });
in
  with import nixpkgs {
    overlays = [(self: super:
    {
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
      cmake
      exa
      gcc
      gdb
      gfortran
      hdf5
      liblapack
      libuuid
      lua5_3
      openmpi
      python3Packages.matplotlib
      python3Packages.numpy
      python3Packages.pandas
      python3Packages.pyyaml
      python3Packages.virtualenvwrapper
      python3Packages.scipy
      valgrind
    ];
  }
