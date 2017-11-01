echo "trees differ across genome with recombination"
python ./examples.py -d 111 -T 1 -N 3 -r 1e-1 -L 10 -a 0.0001 -b 0.0001 -k 3 -t neutral_ts -g /dev/null > /dev/null  && python proc.py
echo " "
echo "trees conserved across genome without recombination"
python ./examples.py -d 111 -T 1 -N 3 -r 0 -L 10 -a 0.0001 -b 0.0001 -k 3 -t neutral_ts -g /dev/null > /dev/null && python proc.py

