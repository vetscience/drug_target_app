sort -t$'\t' -k1,1 -k12,12nr $1 | sort -u -k1,1 --merge > best_hits_$1
