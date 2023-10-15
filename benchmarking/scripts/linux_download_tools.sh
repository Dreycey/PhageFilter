#########################
# Installs for linux VM
#########################

# Basic installs
sudo apt update
sudo apt install libgomp1
sudo apt install zip


# python
sudo apt install python3-pip
sudo apt-get install python3.6
pip install numpy

# JVM for FastViromeExplorer
sudo apt install default-jdk
sudo apt-get install libncurses5-dev libncursesw5-dev
sudo apt install lzma-dev

# for rust cc
sudo apt install build-essential

# install (GLIBCXX)
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install -y g++-11

sudo cp ../linux-binaries/* /bin/
