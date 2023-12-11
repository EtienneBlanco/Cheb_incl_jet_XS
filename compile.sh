g++ Cheb_jet_main.cpp -o jet `lhapdf-config --cflags --ldflags` `TMDlib-config --cppflags --ldflags`
echo "Compiled"
echo
echo "To run the program, use : ./jet"
echo "Parameters can be changed in the file 'Parameters'"
