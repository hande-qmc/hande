let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-unstable";
    sha256 = "1i3p5m0pnn86lzni5y1win0sacckw3wlg9kqaw15nszhykgz22zq";
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
      python36Packages.matplotlib
      python36Packages.numpy
      python36Packages.pandas
      python36Packages.pyyaml
      python36Packages.virtualenvwrapper
      valgrind
    ];
  }
