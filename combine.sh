for f in result*.txt
do
	echo $f >> final.txt
	cat $f >> final.txt
	rm $f
done
