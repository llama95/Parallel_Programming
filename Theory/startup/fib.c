#include <stdio.h>
#include <stdlib.h>


int fibonacci(int term);
int main(){
    int terms = 5;
    int counter;
    for(counter = 0; counter < terms; counter++){
        printf("%d ", fibonacci(counter));
    }

    return 0;
}
// n , n-1 , n-2, n-3, n-4
int fibonacci(int term){
    /* cond to exit recursion*/
    if(term < 2)
       return term;
    return fibonacci(term - 1) + fibonacci(term - 2);
}
