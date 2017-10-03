# include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct Linkedrecords {
  int year;
  int stops;
  struct Linkedrecords *next;
}; struct Linkedrecords

Linkedrecords *list_find(Linkedrecords *head, int year){
  if (head->year == 0){ //base case 1
    head->year = year;
    head->stops = 0;
    head->next = malloc(sizeof(Linkedrecords)); 
    head->next->year = 0;
    head->next->stops = 0;
    head->next->next = 0;

    return head;
  }
  if (head->year == year) //base case 2
    return head;
}
  return(list_find(head->next,year));

int parse_line(char *line){
  char *line_reads = (char *)malloc(sizeof(char)*10);
  strncpy(line_reads,line+6,4);
  int convertToInt = atoi(line_reads); // string to int = atoi
  if (convertToInt){
    return convertToInt;
  }
  return -1;
}
void print_list(Linkedrecords *head){
  if (head == NULL){
    return;
  }
  printf("%d Stops \n",head->year,head->stops);
  print_list(head->next);
}


int main(int argc, char** argv) { // argc number of char* string in our array
  if (argc < 2) {
    printf("error");
}
Linkedrecords *out_list = malloc(sizeof(Linkedrecords)); //allocates size of excel file via out put pointer
FILE *fp;
char *line;
fp = fopen(agrv[1],"r"); // opens our read only file
size_t len =0;
size_t reads;
getline(&line,&len,fp); // no header line needed
Linkedrecords *node = malloc(sizeof(Linkedrecords));
while ((reads = getline(&line,&len,fp) != -1)){
  int year = parse_line(line);
  node = list_find(out_list,year);
  node->stops+=1;

}
fclose(fp);
print_list(out_list);
return 0;



}
