cat $1 | fasta_formatter | awk 'BEGIN{getline; printf $0"\t"; seq="";}{if($0 ~ "^>"){print seq; seq=""; printf $0"\t";}else{seq=seq""$0}}END{print seq}'
