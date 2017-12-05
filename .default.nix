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
      python35Packages.matplotlib
      python35Packages.numpy
      python35Packages.pandas
      python35Packages.virtualenvwrapper
      valgrind
    ];
  }
