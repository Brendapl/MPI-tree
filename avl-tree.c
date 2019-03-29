#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>
typedef enum{FALSE, TRUE} Bool;

/*Node struct with its members.*/
struct node
{
	int val;
	int balance;
	struct node *left_child;
	struct node *right_child;
};

/*Function which searches in a binary tree for a value (data), given a root Node.*/
//   Parallel version   //
/*
struct node* search(struct node *root, int data){
	struct node *foundNode = NULL;
	if(root != NULL){
		if(data == root -> val){
			foundNode = root;
		}
		else if(data < root -> val){
			#pragma omp task shared(foundNode)
			{
				struct node *leftChildNode;
				leftChildNode = search(root -> left_child, data);
				if(leftChildNode != NULL){
					#pragma omp atomic write
						foundNode = leftChildNode;
					#pragma omp cancel taskgroup
				}// end if leftChildNode found
			}//end pragma
		}//end if for left side of subtree
		else if(data > root -> val){
			#pragma omp task shared(foundNode)
			{
				struct node *rightChildNode;
				rightChildNode = search(root -> right_child, data);
				if(rightChildNode != NULL){
					#pragma omp atomic write
						foundNode = rightChildNode;
					#pragma omp cancel taskgroup
				}// end if rightChildNode found
			}//end pragma
		}//end if for right side of subtree
	}// end if root is not NULL
	return(foundNode);
}//end parallel search
*/

/*Function which searches in a binary tree for a value (data), given a root Node.*/
//   Sequential version   //
struct node* search(struct node *root, int data){
	struct node *currNode = root;
	if(currNode != NULL){
		if(data < currNode -> val){
			currNode = search(currNode -> left_child, data);
		}
		else if(data > currNode -> val){
			currNode = search(currNode -> right_child, data);
		}
		else if(currNode -> val != data){
			currNode = NULL;
		}
	}
	return(currNode);
}


struct node* create(int quantity){
	struct node *insert (int data, struct node *ptr, int *ht_inc);
	int rand_max = 10000000;
	int ht_inc;
	int data = rand() % rand_max;

	struct node *root = (struct node *) malloc(sizeof(struct node));

	root -> val = data;
	root -> left_child = NULL;
	root -> right_child = NULL;
	root -> balance = 0;

	for(int i = 1; i < quantity; i++){
		data = rand() % rand_max;
		while(search(root, data) != NULL){
			data = rand() % rand_max;
		}
		root = insert(data, root, &ht_inc);
	}
	return root;
}

double searchFound(struct node *root, int data){
	double elapsed = -1.0;
    struct node *someNode = (struct node *) malloc(sizeof(struct node));

    clock_t start = clock();
    someNode = search(root, data);
    clock_t stop = clock();

		if(someNode != NULL){
    	if(someNode->val == data){

        	elapsed = (double)(stop - start) / CLOCKS_PER_SEC;
        	printf("Time elapsed in ms: %f \n", elapsed);
    	}
	}

    return elapsed;
}

/*Function which inserts in a binary tree and Balance.*/
struct node *insert (int data, struct node *ptr, int *ht_inc)
{
    struct node *aptr;
    struct node *bptr;
    if(ptr==NULL)
    {
        ptr = (struct node *) malloc(sizeof(struct node));
        ptr -> val = data;
        ptr -> left_child = NULL;
        ptr -> right_child = NULL;
        ptr -> balance = 0;
        *ht_inc = TRUE; return (ptr);
    }
    if(data < ptr -> val)
    {
        {
            ptr -> left_child = insert(data, ptr -> left_child, ht_inc);
            if(*ht_inc==TRUE)
                switch(ptr -> balance) {
                    case -1:
                        ptr -> balance = 0;
                        *ht_inc = FALSE;
                        break;
                    case 0:
                        ptr -> balance = 1;
                        break;
                    case 1:
                        aptr = ptr -> left_child;
                        if(aptr -> balance == 1) {
                            printf("Left to Left Rotation\n");
                            ptr -> left_child= aptr -> right_child;
                            aptr -> right_child = ptr;
                            ptr -> balance = 0;
                            aptr -> balance=0;
                            ptr = aptr;
                        }
                        else {
                            printf("Left to Right Rotation\n");
                            bptr = aptr -> right_child;
                            aptr -> right_child = bptr -> left_child;
                            bptr -> left_child = aptr;
                            ptr -> left_child = bptr -> right_child;
                            bptr -> right_child = ptr;
                            if(bptr -> balance == 1 )
                                ptr -> balance = -1;
                            else
                                ptr -> balance = 0;
                            if(bptr -> balance == -1)
                                aptr -> balance = 1;
                            else
                                aptr -> balance = 0;
                            bptr -> balance=0;
                            ptr = bptr;

                        }
                        *ht_inc = FALSE;
                }
        }
    }
    if(data > ptr -> val){
        ptr -> right_child = insert(data, ptr -> right_child, ht_inc);
        if(*ht_inc==TRUE){
            switch(ptr -> balance){
                case 1:
                    ptr -> balance = 0;
                    *ht_inc = FALSE;
                    break;
                case 0:
                    ptr -> balance = -1;
                    break;
                case -1:
                    aptr = ptr -> right_child;
                    if(aptr -> balance == -1) {
                        printf("Right to Right Rotation\n");
                        ptr -> right_child= aptr -> left_child;
                        aptr -> left_child = ptr;
                        ptr -> balance = 0;
                        aptr -> balance=0;
                        ptr = aptr;
                    }
                    else{
                        printf("Right to Left Rotation\n");
                        bptr = aptr -> left_child;
                        aptr -> left_child = bptr -> right_child;
                        bptr -> right_child = aptr;
                        ptr -> right_child = bptr -> left_child;
                        bptr -> left_child = ptr;
                        if(bptr -> balance == -1)
                            ptr -> balance = 1;
                        else
                                ptr -> balance = 0;
                        if(bptr -> balance == 1)
                                    aptr -> balance = -1;
                                else
                                    aptr -> balance = 0;
                        bptr -> balance=0;
                        ptr = bptr;
                    }
                    *ht_inc = FALSE;
            }
        }
    }
    return(ptr);
}

void display(struct node *ptr, int level){
    int i;
    if ( ptr!=NULL ){
        display(ptr -> right_child, level+1);
        printf("\n");
        for (i = 0; i < level; i++)
            printf(" ");
        printf("%d", ptr -> val);
        display(ptr -> left_child, level+1);
    }
}
void inorder(struct node *ptr){
    if(ptr!=NULL){
        inorder(ptr -> left_child);
        printf("%d ",ptr -> val);
        inorder(ptr -> right_child);
    }
}

int * randomFill(int size){
	int i;
	int vec[size];

	for (i = 0; i < size; ++i){
		 vec[i] = i + 1;
	}

	/* Shuffle entries */
//srand( time( 0 ) );
for (i = 0; i < size; ++i) {
		int a = rand( ) % size;
		int b = rand( ) % size;
		if (a != b) {
				int tmp = vec[a];
				vec[a] = vec[b];
				vec[b] = tmp;
		}
	}
	/*for (i = 0; i < size; ++i){

		printf( "%d\n", vec[i]);
	}*/
	return vec;
}

struct node* createTree(int vec[]){
	struct node *insert (int data, struct node *ptr, int *ht_inc);
	int ht_inc;
	struct node *root = (struct node *) malloc(sizeof(struct node));
	root = NULL;
	for(int i = 0; i < 1000000; i++){
		root = insert(vec[i], root, &ht_inc);
	}
	return root;
}


int main(int argc, char *argv[]){

    int ht_inc;
    int data ;
    int option;
		int quantity = 0;
		int valorBuscar = 359688;
		double elapsed;
		double tiempoFunc;
    struct node *root = (struct node *)malloc(sizeof(struct node));
    root = NULL;

		int MASTERID = 0;
		int numNodes = 1000000;

//MPI Variables
		int numCores;
		int resourceCoreId;
		int coreIdDestination;
		int * nodes;
		MPI_Status status;


	//Random Fill Turn to Function
		nodes = randomFill(numNodes);
/*
		for (int i = 0; i < numNodes; ++i){

			printf( "%d\n", nodes [i]);
		}
*/
		coreIdDestination = 1; //CHANGE

		MPI_Init (&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &resourceCoreId);
		MPI_Comm_size(MPI_COMM_WORLD, &numCores);

		int increment = (int) round((double) numNodes / numCores - 1);
		int * start = nodes;
		int finalIncrement = numNodes - increment * (numCores - 2);

		if(resourceCoreId == MASTERID){

			for(int i = 1; i < numCores; i++){
				if(i == numCores - 1){
					//MPI_Send(start, finalIncrement, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
					//MPI_Send(&finalIncrement, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
				} else {
					MPI_Send(start, increment, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
					MPI_Send(&increment, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
					start = (start + increment);
				}
			}
			MPI_Bcast(&valorBuscar, 1, MPI_INT, MASTERID, MPI_COMM_WORLD);

			//MPI_Recv();


		}

		if(resourceCoreId != MASTERID && resourceCoreId != numCores - 1){

			MPI_Recv(&start, increment, MPI_INT, MASTERID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&increment, 1, MPI_INT, MASTERID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			printf("Increment:  %d", increment);
			for (int i = 0; i < increment; ++i){

				printf( "%d\n", start[i]);
			}

			/*
			root = createTree(nodes);

			tiempoFunc = searchFound(root, valorBuscar);

			MPI_Send(&tiempoFunc, 1, MPI_DOUBLE, MASTERID, MPI_ANY_TAG, MPI_COMM_WORLD);
			*/
		}
		MPI_Finalize();

		//int time = searchFound(elapsed);
}
