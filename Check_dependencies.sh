#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Function description
# Check dependencies for 3dToFunc 
#===============

# check the R version
Rversion=$(R --version | head -n 1 | awk -F[" "] '{print $3}' -)
if [[ -z "$Rversion" ]]; then
    echo "ERROR ====>>> No R installation is detected in the system !!!" 
    errcond=1
fi
numfield=`echo $Rversion | awk -F'[.]' '{print NF}' -`
if [[ $numfield -ge 3 ]]; then
	num1=`echo $Rversion | awk -F'[.]' '{print $1}' -`
	if [[ $num1 -ge 4 ]]; then
		echo "Installed R version: "${Rversion}
	else 
		echo " -- R version should be at least R 4.0 !!! "
		errcond=1
	fi
else
	num1=`echo $Rversion | awk -F'[.]' '{print $1}' -`
	num2=`echo $Rversion | awk -F'[.]' '{print $2}' -`
	if [[ $num1 -ge 4 ]]; then
		echo "Installed R version: "${Rversion}
	else 
		echo " -- R version should be at least R 4.0 !!!"
		errcond=1
	fi	
fi

# check the g++ version
if command -v g++ >/dev/null 2>&1 ; then
    echo "Installed g++"
    cd ./source/C++/
    g++ -std=c++0x -o straw main.cpp straw.cpp -lcurl -lz
    if [ -f "./straw" ]; then
        echo "Checking dependencies successfully!"
    fi
else
    echo "g++ not found"
fi
