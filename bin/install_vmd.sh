#!/bin/bash
#0. tested on mac os 10.13.6, 10.14
#1. install Xcode first
#2. Download vmd source v 1.9.3
#3. Add vmdmac.r to the src directory https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/0678.html
#Set OSX version manually! - that's important for fltk
export OSX_VER=10.13.6
#compiling fltk
#We need latest development version
cd vmd-1.9.3/lib/fltk
wget http://fltk.org/pub/fltk/snapshots/fltk-1.4.x-r13117.tar.gz
tar -xf fltk-1.4.x-r13117.tar.gz
ln -s fltk-1.4.x-r13117 fltk
ln -s fltk include
cd fltk
./configure --prefix="$PWD/../MACOSXX86_64" --exec-prefix="$PWD/../MACOSXX86_64" --libdir="$PWD/../MACOSXX86_64" CXXFLAGS="-mmacosx-version-min=$OSX_VER" LDFLAGS="-mmacosx-version-min=$OSX_VER"
make -j 8
make install
cd ..
cd ..
cd ..
cd ..
export PLUGINDIR="$PWD/vmd-1.9.3/plugins"
export export TCLINC=-I/System/Library/Frameworks/Tcl.framework/Versions/8.5/Headers
export export TCLLIB=-L/System/Library/Frameworks/Tcl.framework/Versions/8.5/Headers
cd plugins
make   MACOSXX86_64 TCLINC=$TCLINC TCLLIB=$TCLLIB
make   distrib 
cd ../vmd-1.9.3
echo "MACOSXX86_64 LP64 FLTKOPENGL FLTK TK  TCL PTHREADS " > configure.options
#Fix code
sed -i.bak 's/MACOSX/MACOSXX86/g' bin/vmd.sh
sed -i.bak 's/MACOSX/MACOSXX86/g' bin/vmd.csh
export VMDINSTALLBINDIR=$PWD/../../vmd-bin/bin #/usr/local/bin
export VMDINSTALLLIBRARYDIR=$PWD/../../vmd-bin/vmd #/usr/local/lib/$install_name
./configure 
cd src
sed -i.bak 's/fltk-1.3.x/fltk/g' Makefile
sed -i.bak 's%../lib/tk/lib_MACOSXX86_64/Tk.framework/Versions/8.5/Headers%/System/Library/Frameworks/Tk.framework/Versions/8.5/Headers%g' Makefile
make veryclean
make -j 8
make install
$VMDINSTALLBINDIR/vmd

