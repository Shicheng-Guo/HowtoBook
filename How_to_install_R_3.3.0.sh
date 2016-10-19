

./configure --prefix=$HOME/software/R-3.3.0 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' \
 '--with-blas' '--with-lapack' '--enable-R-profiling' \
 '--enable-R-shlib' \
 '--enable-memory-profiling' '--enable-R-static-lib'
 make clean 
 make
 
 
 
  cd ~/software
  wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
  tar xzvf bzip2-1.0.6.tar.gz
  cd bzip2-1.0.6
  make -n install PREFIX=$HOME/software/bzip2-1.0.6
  make install PREFIX=$HOME/software/bzip2-1.0.6
  make -f Makefile-libbz2_so
  
  ln -s libbz2.so.1.0.6 libbz2.so.1.0

  
 ../configure --prefix=$HOME/software/R-3.3.0 '--with-cairo' \
 '--with-jpeglib' '--with-readline' '--with-tcltk' \
 '--with-blas' '--with-lapack' '--enable-R-profiling' \
 '--enable-R-shlib' \
 '--enable-memory-profiling' \  '--enable-R-static-lib'
 
 
 cd ~/software
 wget http://tukaani.org/xz/xz-5.2.2.tar.gz
 tar xzvf xz-5.2.2.tar.gz
 cd xz-5.2.2
 ./configure --prefix=$HOME/software/xz-5.2.2	
 make -j3
 make install

 # symlink libbz2.so.1 to libbz2.so.1.0.6
  ln -s libbz2.so.1.0.6 libbz2.so.1.0
  ln -s libbz2.so.1.0 libbz2.so.1
  ln -s libbz2.so.1 libbz2.so

libbz2.so.1.0.6 to libbz2.so.1.0
libbz2.so.1.0 to libbz2.so.1
libbz2.so.1 to libbz2.so



 cd ~/software
 wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.38.tar.gz
 tar xzvf pcre-8.38.tar.gz
 ./configure --prefix=$HOME/software/pcre-8.38
 make -j3
 make install
 
 
 cd ~/software
 wget ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.38.tar.gz
 
 tar xzvf pcre-8.38.tar.gz
 cd pcre-8.38
 ./configure --enable-utf8 --prefix=$HOME/software/pcre-8.38
 make -j3
 make install

 
 cd ~/software
 wget --no-check-certificate https://curl.haxx.se/download/curl-7.47.1.tar.gz
 tar xzvf curl-7.47.1.tar.gz
 cd curl-7.47.1
 ./configure --prefix=$HOME/software/curl-7.47.1
 make -j3
 make install
 

 
 
export PATH=/home/shg047/software/pcre-8.38/bin:$PATH
export PATH=/home/shg047/software/xz-5.2.2/bin:$PATH
export PATH=/home/shg047/software/zlib-1.2.8:$PATH
export PATH=/home/shg047/software/bzip2-1.0.6:$PATH
export PATH=/home/shg047/software/curl-7.47.1/bin:$PATH

export CFLAGS="-I/home/shg047/software/bzip2-1.0.6 -I/home/shg047/software/zlib-1.2.8 -I/home/shg047/software/xz-5.2.2/include -I/home/shg047/software/pcre-8.38/include -I/home/shg047/software/curl-7.47.1/include"
export LDFLAGS="-L/home/shg047/software/bzip2-1.0.6 -L/home/shg047/software/zlib-1.2.8 -L/home/shg047/software/xz-5.2.2/lib -L/home/shg047/software/pcre-8.38/lib -L/home/shg047/software/curl-7.47.1/lib"

export CFLAGS="-I/home/shg047/software/xz-5.2.2/include"
export LDFLAGS="-L/home/shg047/software/xz-5.2.2/lib"



./configure --enable-utf8 --prefix=/home/shg047/software/


