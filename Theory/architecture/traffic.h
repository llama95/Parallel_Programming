#ifndef traffic_h
#define traffic_h

typedef struct Linkedrecords {
  int year;
  int stops;
  struct Linkedrecords *next;
} Linkedrecords;

Linkedrecords *list_find(Linkedrecords *head, int year);
int parse_line(char *line);
void print_list(Linkedrecords *head);
#endif /* traffic_h */
