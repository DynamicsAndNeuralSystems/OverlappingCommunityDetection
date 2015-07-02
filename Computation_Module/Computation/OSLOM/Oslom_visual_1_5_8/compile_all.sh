

echo "Compiling OSLOm undirected (oslom_undir) ..."
make -f Sources/Makefiles/makefile_undir

echo ""
echo "Compiling OSLOm directed (oslom_dir) ..."
make -f Sources/Makefiles/makefile_dir



echo ""
echo "Compiling program to write pajek format (pajek_write_undir) ..."
make -f Sources/Makefiles/makefile_pajek_undir


echo ""
echo "Compiling program to write pajek format (pajek_write_dir) ..."
make -f Sources/Makefiles/makefile_pajek_dir


echo ""
echo "Compiling infomap_undirected ..."
cd Sources/infomap_undir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv Sources/infomap_undir/infomap infomap_undir
mv Sources/infomap_undir/infomap_scr infomap_undir_script


echo ""
echo "Compiling infomap_directed ..."
cd Sources/infomap_dir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv Sources/infomap_dir/infomap infomap_dir
mv Sources/infomap_dir/infomap_scr infomap_dir_script


echo ""
echo "Compiling louvain  method ..."
cd Sources/louvain/
g++ script_to_compile.cpp
./a.out
g++ -o louvain_script order.cpp -O3
cd ../..
mv Sources/louvain/louvain_script .
mv Sources/louvain/community .
mv Sources/louvain/convert .
mv Sources/louvain/hierarchy .





echo ""
echo "***********************************************"
echo "type \"./oslom_undir\" (or \"./oslom_dir\") to check if the program is compiled and works"
