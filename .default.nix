let
  hostPkgs = import <nixpkgs> {};
  nixpkgs = (hostPkgs.fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "nixos-18.03";
    sha256 = "1g4vfvpvdbgap84qkngzcai2gqwbpv77sqyswv884k1qxf2xji7j";
  });
in
  with import nixpkgs {
    overlays = [(self: super:
      {
        hdf5 = super.hdf5.override {
          gfortran = super.gfortran;
          mpi = super.openmpi;
        };
        python3 = super.python3.override {
          packageOverrides = py-self: py-super: {
            matplotlib = py-super.matplotlib.override {
              enableTk = true;
              enableQt = true;
            };
          };
        };
      }
    )];
  };
  stdenv.mkDerivation {
    name = "HANDE";
    buildInputs = [
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
      pipenv
      python3Full
      python3Packages.matplotlib
      python3Packages.numpy
      python3Packages.pandas
      python3Packages.pyyaml
      python3Packages.scipy
      stdenv
      valgrind
    ];
    src = null;
    shellHook = ''
    SOURCE_DATE_EPOCH=$(date +%s)
    '';
  }
