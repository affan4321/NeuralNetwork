#include <iostream>
#include <string>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <fstream>
#include <cstring>
#include <sstream>

using namespace std;

struct nodeData{
    int id;
    float** weights;
    float finalWeight;
    float* I;
};

sem_t s1;

string _pipe2 = "pipe2"; 
string _pipe3 = "pipe3";

void* hiddenLayer(void* args){
    nodeData* l = (nodeData*)args;

    for(int i = 0; i < 8; i++){
        if(l->id == i+1){
            l->finalWeight = (l->I[0] * l->weights[0][i]) + (l->I[1] * l->weights[1][i]);
        }
    }

    sem_post(&s1);
    // writing to pipe 3
    cout << "writing pipe 3" << endl;
    int fd;
    fd = open(_pipe3.c_str(), O_WRONLY);
    write(fd, l, sizeof(l));
    close(fd);

    pthread_exit(NULL);
}


int main(int argc, char* argv[])
{
    sem_init(&s1, 0, 0);
    cout << "Entered hidden layer" << ++*(argv[1]) << endl;
    nodeData l;
    cout << "reading pipe 2" << endl;
    // Reading data from pipe 2
    int fd;
    fd = open(_pipe2.c_str(),O_RDONLY);
    read(fd, &l, sizeof(l));
    close(fd);

    // Making threads of size of weights provided;
    pthread_t* h = new pthread_t[8];
    // set thread attributes
    pthread_attr_t* attr = new pthread_attr_t[8];

    cout << "creating threads" << endl;
    for(int i = 0; i < 8; i++){
        pthread_attr_init(&attr[i]);
        pthread_attr_setdetachstate(&attr[i], PTHREAD_CREATE_DETACHED);
        // set thread scheduling policy
        pthread_attr_setschedpolicy(&attr[i], SCHED_FIFO);
        // create thread
        // assigning id to thread
        l.id = i+1;
        pthread_create(&h[i], &attr[i], hiddenLayer, &l);
    }

    //sem_wait(&s1);

    // reading from pipe 3
    cout << "reading pipe 3" << endl;
    nodeData n;
    fd = open(_pipe3.c_str(), O_RDONLY);
    read(fd, &n, sizeof(n));
    close(fd);

    cout << "node data read: " << endl;
    for(int i = 0; i < 8; i++){
        cout << n.finalWeight << " ";
    }
    cout << endl;

    cout << "Program reached at the end of hidden Layer" << endl;

    return *argv[1];
}


/*

Creation of a neural network with 3 layers:

Input layer:
input  -> 2 weights
arr[2] -> weight1, weight2
pipe1  -> arr[2]
call hidden layer

Hidden layer:
    ForwardLayer
            -->Forward propogation
    read from pipe1
    create Matrix
    pipe2( use between ForwardLayer and CalcculationLayer)

    CalculationLayer:
            -->Backpropagation:
    {0,1,2,3,4}-> exec -> process -> mat[0]->thread mat[1]->thread...
    (communication between files through pipes)

    pipe3(use between CalcculationLayer and ReturnResultLayer)
    ReturnResultLayer:
        calculated gradients -> f(x) -> send to output layer

    pipe4(use between ReturnResultLayer and OutputLayer)

Output layer:
        weights -> output

*/


/*
    1- Neuron Structure:

struct Neuron {
    double input;
    double output;
    std::vector<double> weights;
}; 
2- Layer Structure:

struct Layer {
    std::vector<Neuron> neurons;

    void initialize(int numNeurons, int numInputs) {
        neurons.resize(numNeurons);
        for (auto& neuron : neurons) {
            neuron.weights.resize(numInputs);
        }
    }

    void forwardPropagate(const std::vector<double>& inputs, std::vector<double>& outputs) {
        outputs.resize(neurons.size());
        for (size_t i = 0; i < neurons.size(); ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < inputs.size(); ++j) {
                sum += neurons[i].weights[j] * inputs[j];
            }
            outputs[i] = sum;
        }
    }
}; 
3 - Neural Network Initialization:


struct NeuralNetwork {
    std::vector<Layer> layers;

    void addLayer(int numNeurons, int numInputs) {
        Layer layer;
        layer.initialize(numNeurons, numInputs);
        layers.push_back(layer);
    }

    void initialize() {
        for (size_t i = 1; i < layers.size(); ++i) {
            int numInputs = layers[i - 1].neurons.size();
            layers[i].initialize(layers[i].neurons.size(), numInputs);
        }
    }

    void propagateInputs(const std::vector<double>& inputs, std::vector<double>& outputs) {
        for (size_t i = 0; i < layers.size(); ++i) {
            if (i == 0) {
                layers[i].forwardPropagate(inputs, outputs);
            } else {
                layers[i].forwardPropagate(outputs, outputs);
            }
        }
    }
}; 
4 - Multi-threading and Communication:


// Example code for two layers running in separate processes
// Layer 1 Process
void layer1Process(int readPipe, int writePipe) {
    // Read inputs from the pipe
    std::vector<double> inputs;
    // Read inputs from the readPipe

    // Perform computations
    std::vector<double> outputs;
    neuralNetwork.layers[0].forwardPropagate(inputs, outputs);

    // Send outputs to the next layer through the writePipe
    // Write outputs to the writePipe
}

// Layer 2 Process
void layer2Process(int readPipe, int writePipe) {
    // Read inputs from the pipe
    std::vector<double> inputs;
    // Read inputs from the readPipe

    // Perform computations
    std::vector<double> outputs;
    neuralNetwork.layers[1].forwardPropagate(inputs, outputs);

    // Send outputs to the next layer through the writePipe
    // Write outputs to the writePipe
} 
 5 - Forward Propagation:


std::vector<double> inputs = {1.0, 2.0, 3.0}; // Example inputs
std::vector<double> outputs;

// Propagate inputs through the layers
neuralNetwork.propagateInputs(inputs, outputs); 
6 - Backward Propagation:


// Example backward propagation from the output layer to the first layer
for (size_t i = neuralNetwork.layers.size() - 1; i > 0; --i) {
    // Read inputs from the pipe
    std::vector<double> inputs;
    // Read inputs from the readPipe

    // Send outputs to the previous layer through the
#include <unistd.h> // For pipe, fork, close
#include <sys/wait.h> // For wait

void layer1Process(int readPipe, int writePipe) {
    // Close unused pipe ends
    close(readPipe);

    // Read inputs from the pipe
    std::vector<double> inputs;
    // Read inputs from the writePipe

    // Perform computations
    std::vector<double> outputs;
    neuralNetwork.layers[0].forwardPropagate(inputs, outputs);

    // Send outputs to the next layer through the writePipe
    // Write outputs to the readPipe

    // Close the pipe after writing
    close(writePipe);
}
void layer2Process(int readPipe, int writePipe) {
    // Close unused pipe ends
    close(writePipe);

    // Read inputs from the pipe
    std::vector<double> inputs;
    // Read inputs from the readPipe

    // Perform computations
    std::vector<double> outputs;
    neuralNetwork.layers[1].forwardPropagate(inputs, outputs);

    // Send outputs to the next layer through the writePipe
    // Write outputs to the writePipe

    // Close the pipe after writing
    close(readPipe);
}

int main() {
    int pipe1[2]; // Pipe between layer 1 and layer 2
    int pipe2[2]; // Pipe between layer 2 and layer 1

    // Create the pipes
    if (pipe(pipe1) == -1 || pipe(pipe2) == -1) {
        perror("Pipe creation failed");
        exit(1);
    }

    pid_t child1 = fork();
    if (child1 == -1) {
        perror("Fork failed");
        exit(1);
    } else if (child1 == 0) {
        // Child process (layer 1)
        layer1Process(pipe1[0], pipe2[1]);
        exit(0);
    }

    pid_t child2 = fork();
    if (child2 == -1) {
        perror("Fork failed");
        exit(1);
    } else if (child2 == 0) {
        // Child process (layer 2)
        layer2Process(pipe2[0], pipe1[1]);
        exit(0);
    }

    // Parent process

    // Close unused pipe ends
    close(pipe1[0]);
    close(pipe1[1]);
    close(pipe2[0]);
    close(pipe2[1]);

    // Wait for child processes to finish
    waitpid(child1, nullptr, 0);
    waitpid(child2, nullptr, 0);

    return 0;
}
*/
