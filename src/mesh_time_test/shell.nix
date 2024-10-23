{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    (pkgs.python3.withPackages (ps: with ps; [
      statistics
      numpy
      matplotlib
    ]))
  ];

  shellHook = ''
    alias \
      q="exit" \
      cl="clear" \
      nv="nvim"
  '';
}
