awk 'BEGIN{FS=",";} {
  for(i = 1; i <= NF; i++) {
    if(match($i, /[A-Za-z]/)) printf("%s ", $i);
    else printf("%.2f ", $i);
    if(i < NF) printf("& ");
    else printf("\\\\ \n");
  }
}' $1
