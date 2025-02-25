make
if [ "$?" != 0 ]
then
	echo "make error"
    exit 1 #参数正确，退出状态为0
fi
echo "start "
ulimit -s unlimited
ulimit -v unlimited
cat small.gra | ./main


