Optional: Permanently enable scl toolchain by putting this in your .bashrc (warning: don't try to use the scl enable devtoolset-2 bash command from before in your .bashrc. This spawns a new bash shell, and if that's in your .bashrc, it creates a new shell, which loads your .bashrc, which creates a new shell, etc.)
source /opt/rh/devtoolset-2/enable
may have worked also
