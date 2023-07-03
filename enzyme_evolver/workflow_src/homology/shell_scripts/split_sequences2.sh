cat $1 |sed 's/ //g;s/\[//g;s/\]//g;s/\;//g;s/\;//g;s/\|//g;/^$/d'| awk -v dir=$2 '{
if (substr($0, 1, 1)==">") {filename=dir"/"(substr($0,2)) ".fasta"; print(substr($0,2))}
                        print $0 > filename
                }'
