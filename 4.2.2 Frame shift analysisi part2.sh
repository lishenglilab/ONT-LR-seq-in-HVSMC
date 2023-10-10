work_dir='/home/haowu/project/pinghuaji/full_lengh_pinghuaji/Result/ORFfind'
cd $work_dir
conditions=`ls | grep condition`
for condition in $conditions
do
	samples=`ls -lh $condition | grep d | awk '{print$9}'`
	for i in $samples
	do
		makeblastdb -in $condition/$i/ref.fa -dbtype prot -out $condition/$i/database
		blastp -db $condition/$i/database -query $condition/$i/query.fa -outfmt 6 -out $condition/$i/query.res.txt -num_threads 10
		sed -i '1 i\qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' $condition/$i/query.res.txt
	done
done

for condition in $conditions
do
	samples=`ls -lh $condition | grep d | awk '{print$9}'`
	for i in $samples
	do
		blastx -db $condition/$i/database -query $condition/$i/query.cds.fa -outfmt 0 -out $condition/$i/query.res.frameshift.txt -num_threads 10
	done
done