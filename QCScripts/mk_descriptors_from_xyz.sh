#!/bin/tcsh

set NPROC = `nproc`

set O_WORKDIR = `pwd`

mkdir -p /dir/to/QCSCR
setenv QCSCRATCH /dir/to/QCSCR

alias print_time 'set t=`date` && echo "[$t] $1"'

echo "═══════════════════════════════════════════════"
print_time "START"

# 1
set Scriptdir = "/dir/to/QCScripts/"
cp $Scriptdir"/mk_opt_in.py" "./"
cp $Scriptdir"/opt_format.in" "./"
python mk_opt_in.py
rm mk_opt_in.py
rm opt_format.in

# 2
cp $Scriptdir"/mk_rst_in.py" "./"
cp $Scriptdir"/rst_format.in" "./"
cp $Scriptdir"/read_rst_data.py" "./"

# 3
if (! -e "descriptor.csv") then
	touch "descriptor.csv"
endif

# iter with xyz files
set target_text = "PCnamehere"
cd $O_WORKDIR/xyz
foreach file (*.xyz)
	set PCname = `basename $file .xyz`
	echo "═══════════════════════════════════════════════"
	echo "Processing file: $PCname"
	
	# 4
	echo "----------------------------------------------------------------"
	print_time "Geometry optimization start"
	cd $O_WORKDIR/opt
	set optname = $PCname"_opt"
	if(!(-e $optname.out) && (-e $optname.in)) then
		qchem -np $NPROC $optname.in $optname.out
		if ($status == 0) then
			rm $optname.in
			print_time "Geometry optimizatoin done"
		else
			print_time "Error: Geometry optimization did not finish successfully for $PCname"
		endif
	else
		if((-e $optname.out) && (-e $optname.in)) then
			rm $optname.in
			print_time "Geometry optimization had already been completed"
		endif
	endif

	
	# 5
	echo "----------------------------------------------------------------"
	print_time "reorgE, dE, T1, S1 calculation start"
	cd $O_WORKDIR
	sed -i "s/$target_text/$PCname/g" "./mk_rst_in.py"
	python mk_rst_in.py
	sed -i "s/$PCname/$target_text/g" "./mk_rst_in.py"
	cd $O_WORKDIR/rst
	set rstname = $PCname"_rst"
	if(!(-e $rstname.out) && (-e $rstname.in)) then
		qchem -np $NPROC $rstname.in $rstname.out
		if ($status == 0) then
			rm $rstname.in
			print_time "reorgE, dE, T1 calculation done"
		else
			print_time "Error: reorgE, dE, T1 calculation did not finish successfully for $PCname"
		endif
	endif
	echo "----------------------------------------------------------------"

	# 6
	cd $O_WORKDIR
	sed -i "s/$target_text/$PCname/g" "./read_rst_data.py"
	python read_rst_data.py
	sed -i "s/$PCname/$target_text/g" "./read_rst_data.py"
end

# 7
cd $O_WORKDIR
rm mk_rst_in.py
rm rst_format.in
rm read_rst_data.py

echo "═══════════════════════════════════════════════"
print_time "FINISH"
echo "═══════════════════════════════════════════════"
