set -e
S=../src; M="python3 memtime.py"
python3 gen.py big 5 110 2000000 3
$M $S/rz-prep big/list.tsv big
$M /root/lmo/deps/Big-BWT/bigbwt big.S > big.bigbwt.log
$M $S/rz-lz77 -t 2 -w 64M -b 64M big.S big.left
$M $S/rz-lz77 -r -t 2 -w 64M -b 64M big.S big.right
$M $S/rz-build -a big.S big.S.bwt big.left big.right big.tbl big
for m in 10 20 50 100 200; do $S/rz-genpat big.S $m 1000 pbig$m.txt $m; $S/rz-bench big.rz big.rix pbig$m.txt; done
echo ALLDONE
