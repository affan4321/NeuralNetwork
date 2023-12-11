#include "neuralNetwork.h"
#include <semaphore.h>

using namespace std;

NeuralNetwork n;
int threadCount = 0;


sem_t s,s2;

pthread_mutex_t lock;

void * processThread(void * args)
{
    pthread_mutex_lock(&lock);
    cout << "thread called with id: " << ++threadCount << endl;

    NeuralNetwork neural = *(NeuralNetwork*) args;

    int numOfneurons = 0,outputNeurons = 0;
    int currentLayer = neural.currentLayer;
    float* input;
    float* output;

    //defining an input size to read
    if(currentLayer == 0)
    {   numOfneurons = neural.neurons_initial;
        outputNeurons = neural.neurons_hidden;
        output = new float[outputNeurons];
        input = new float [numOfneurons];
    }

    else if(currentLayer>0 && currentLayer != n.numOflayers-1)
    {
        numOfneurons = neural.neurons_hidden;
        input = new float[numOfneurons];
        output = new float[numOfneurons];
    }
    else if(currentLayer == n.numOflayers-1){
        numOfneurons = neural.neurons_hidden;
        input = new float[numOfneurons];
        output = new float;
    }

    // Reading input pipe
    cout << "\nreading inputPipe "<<currentLayer+1<<"\n" << endl;
    read(neural.inputPipe[currentLayer][0],input,sizeof(float)*neural.layers[currentLayer].numOfNeurons);

    cout<<"inputPipe = ";
    for(int i=0;i<numOfneurons;i++)
        cout << input[i]<< ",";
    cout <<endl;

    //input layer calculations
    if(currentLayer == 0)
    {
        // Calculate sum for each neuron in the input layer
        for (int i = 0; i <neural.layers[1].numOfNeurons; i++) 
        {
            float sum = 0.0;

            for (int j = 0; j < neural.layers[0].numOfNeurons; j++) 
            {
                // Multiply the input value with the corresponding weight from the input layer
                sum += input[j] * neural.layers[0].weights[j][i];
            }
            // Append the resulting sum to the outputs vector
            output[i] = sum;
        }
        cout << endl;
    }

    else if(currentLayer>0 &&  currentLayer != n.numOflayers-1)
    {
        // Calculate sum for each neuron in the input layer
        for (int i = 0; i <neural.layers[currentLayer].numOfNeurons; i++) 
        {
            float sum = 0.0;
            for (int j = 0; j < neural.layers[currentLayer].numOfNeurons; j++) 
            {
                // Multiply the input value with the corresponding weight from the input layer
                sum += input[j] * neural.layers[currentLayer].weights[j][i];
            }
            // Append the resulting sum to the outputs vector
            output[i] = sum;
        }
        cout << endl;
    }
    else if(currentLayer == n.numOflayers-1){
        
        float sum = 0.0;
        for (int j = 0; j < neural.layers[currentLayer].numOfNeurons; j++) 
        {
            // Multiply the input value with the corresponding weight from the input layer
            sum += input[j] * neural.layers[currentLayer].weights[0][j];
        }
        // Append the resulting sum to the outputs vector
        output[0] = sum;
    }

    // Showing contents of inputPipe of upcoming layer
    cout<<"Output sent to pipe = ";
    if(currentLayer != n.numOflayers-1){
        for(int i=0;i<neural.neurons_hidden;i++)
        {
            cout << output[i]<< ",";
        }
        cout <<endl;
    }
    else if(currentLayer == n.numOflayers-1){
        cout << output[0];
    }
    cout << endl;
    
    // Updating inputPipe and making route for upcoming Layer
    if (currentLayer + 1 < neural.numOflayers) 
    {
        n.currentLayer++;
        cout << "\nwriting inputPipe "<<n.currentLayer+1<<"\n" << endl;
        write(neural.inputPipe[n.currentLayer][1], output, sizeof(float)*neural.neurons_hidden);

    }
    else{
        float backPropagate[2];
        backPropagate[0]= ((output[0] * output[0]) + output[0]+1)/2;
        backPropagate[1]= ((output[0] * output[0]) - output[0])/2;

        cout << "backpropagation Values: " << backPropagate[0]<< "," << backPropagate[1] <<endl; 

        cout << "writing to output pipe" << endl;
        write(neural.outputPipe[currentLayer][1], backPropagate, sizeof(float)*2);
    }

    cout << "Thread with id: " << threadCount << " Exiting" << endl;
    pthread_mutex_unlock(&lock);
    pthread_exit(NULL);
}


void* backpropagate(void* args) {
    
    pthread_mutex_lock(&lock);

    NeuralNetwork neural = *(NeuralNetwork*) args;
    cout << "\nbackTracking thread id: " << threadCount-- << endl;
    int currentLayer=neural.currentLayer;
    
    float backPropagate[2];

    // Reading the contents of outputPipe
    cout << "reading outputPipe" << endl;
    read(neural.outputPipe[currentLayer][0], backPropagate, sizeof(float)*2);

    cout << "backpropagation Values: " << backPropagate[0]<< "," << backPropagate[1] <<endl; 

    // Writing to outputPipe
    cout << "writing outputPipe" << endl;
    if(currentLayer!=0)
    {
        write(neural.outputPipe[currentLayer-1][1], backPropagate, sizeof(float)*2);
        n.currentLayer--;
    }
    else if(currentLayer==0)
        write(neural.outputPipe[0][1], backPropagate, sizeof(float)*2);

    cout << "Thread with id: " << threadCount+1 << " Exiting" << endl;
    sem_post(&s2);
    pthread_mutex_unlock(&lock);

    pthread_exit(NULL);
} 

int main(){

    int layers;
    cout << "enter how many layers: ";
    cin >> layers;


    float input[2]= {0.1,0.2};


    n.initialize(layers,2,8,8);
    n.readInputs();
    //n.displayfilesData(layers);       // To display that all initial weights are read perfectly
    pthread_mutex_init(&lock,NULL);

    pthread_t pid[layers];      // Process Threads
    pthread_t pid2[layers];     // BackTrack Threads
    float finalOutput[2];

    for(int k=0;k<2;k++)
    {
        sem_init(&s, 0, layers);
        sem_init(&s2, 0, layers);
        cout<<"Pass number " << k+1 << endl << endl;
        cout << "\nwriting inputPipe 1\n" << endl;
        write(n.inputPipe[0][1],&input,sizeof(input));

        // Initializing mutex lock variable

        // Creating Process Threads
        for(int i=0;i<layers;i++)
            pthread_create(&pid[i],NULL,processThread,(void*)&n);
        
        for(int i=0;i<layers;i++)
        {
            pthread_join(pid[i],NULL);
        }

        // Creating BackTrack Threads
        for(int i=0;i<layers;i++)
            pthread_create(&pid2[i],NULL,backpropagate,(void*)&n);

        for(int i=0;i<layers;i++)
            pthread_join(pid2[i],NULL);
        
        // Waiting on Threads to complete before exit
        for(int i = 0; i < n.numOflayers; i++)
            sem_wait(&s);
        
        for(int i = 0; i < n.numOflayers; i++)
            sem_wait(&s2);

        
        cout << "\nreading outputPipe final value" << endl;
        // Reading first value of Output pipe that is backtracked
        read(n.outputPipe[0][0], finalOutput, sizeof(float)*2);

        cout << "finalOutput: " << finalOutput[0]<< "," << finalOutput[1] <<endl;

        input[0]=finalOutput[0];
        input[1]=finalOutput[1];

    }

     pthread_mutex_destroy(&lock);
    
    // Destroy Mutex
    return 0;
}