/*
I tested my code for M/m/1 queue as well and as talked with the professor, if my code works for M/m/1 queue I should be getting marks for HW3 as well
Thanks!!
Sina
*/

/*
One more thing to note I have implemented the logic for start_service in my arrival and departure as I found it easier that way to keep track of the time and queue for the simulated values
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ARRIVAL 1
#define DEPARTURE 3

double exponentialRandom(double rate) {
    return -log(1.0 - ((double)rand() / RAND_MAX)) / rate;
}

typedef struct Event {
    double time;         // Time of the event (arrival or departure)
    int type;            // 1 for Arrival, 3 for Departure
    double arrivalTime;  // Store arrival time for calculating response time at departure
    struct Event* next;
} Event;

typedef struct {
    Event* head;
} EventList;

typedef struct {
    double arrivalRate;
    double serviceRate;
    int servers;
    int maxQueueSize;
    int customerServed;
    int customerDropped;
    double currentTime;
    int totalCustomers;
    double totalServiceTime;
    double totalQueueLength;
    double totalResponseTime;
	int currentQueueLength; // New member for actual queue management
    double totalWaitingTime;
    double idleTime; // Time when no customers are being served
    double lastEventTime;
    int inService; // Number of customers currently being served
	double* waitingTimes; // Array to store arrival times of waiting customers
    int waitingCapacity;  // Capacity of the waitingTimes array
    int waitingCount;     // Number of waiting customers
	double idleStartTime; // Time when the system last became idle
	double computedEn;
    double computedEr;
    double computedEw;
    double computedP0;
    double computedLambdaEff;
} Simulation;

Simulation* initializeSimulation(double lambda, double mu, int m, int B, int D, int seed) {
    Simulation* sim = (Simulation*)malloc(sizeof(Simulation));
    if (!sim) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    sim->arrivalRate = lambda;
    sim->serviceRate = mu;
    sim->servers = m;
    sim->maxQueueSize = B;
    sim->customerServed = 0;
    sim->customerDropped = 0;
    sim->currentTime = 0.0;
    sim->totalCustomers = D;
	sim->currentQueueLength = 0; // Initialize actual queue length
    sim->totalServiceTime = 0.0;
    sim->totalQueueLength = 0.0;
    sim->lastEventTime = 0.0;
	sim->totalWaitingTime = 0.0;
	sim->totalResponseTime = 0.0;
	sim->idleTime = 0.0;
    sim->inService = 0;
	sim->waitingTimes = NULL;
	sim->waitingCapacity = 0;
	sim->waitingCount = 0;
	sim->idleStartTime = -1; // -1 indicates no idle period has started
	
	sim->computedEn = 0.0;
    sim->computedEr = 0.0;
    sim->computedEw = 0.0;
    sim->computedP0 = 0.0;
    sim->computedLambdaEff = 0.0;
	
    srand(seed); // Initialize the random number generator
    return sim;
}

EventList* createEventList() {
    EventList* list = (EventList*)malloc(sizeof(EventList));
    if (!list) {
        fprintf(stderr, "Memory allocation error for EventList\n");
        exit(EXIT_FAILURE);
    }
    list->head = NULL;
    return list;
}

Event* createEvent(double time, int type) {
    Event* newEvent = (Event*)malloc(sizeof(Event));
    if (!newEvent) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    newEvent->time = time;
    newEvent->type = type;
    newEvent->next = NULL;
    return newEvent;
}

void insertEvent(EventList* list, Event* event) {
    if (list->head == NULL || list->head->time > event->time) {
        // Insert at the beginning
        event->next = list->head;
        list->head = event;
    } else {
        // Find the correct position
        Event* current = list->head;
        while (current->next != NULL && current->next->time <= event->time) {
            current = current->next;
        }
        event->next = current->next;
        current->next = event;
    }
}


void processArrival(Simulation* sim, Event* event, EventList* list) {
    double timeElapsed = event->time - sim->lastEventTime;
    sim->totalQueueLength += sim->currentQueueLength * timeElapsed;
    sim->lastEventTime = event->time;

    // Check if the system was idle before this arrival
    if (sim->inService == 0 && sim->idleStartTime >= 0) {
        sim->idleTime += event->time - sim->idleStartTime; // Update total idle time
        sim->idleStartTime = -1; // Reset since system is no longer idle
    }

    sim->currentTime = event->time;
    if (sim->inService < sim->servers) {
        // No need to wait, service starts immediately
        sim->inService++;
        double serviceTime = exponentialRandom(sim->serviceRate);
        Event* departureEvent = createEvent(sim->currentTime + serviceTime, DEPARTURE);
        departureEvent->arrivalTime = event->time; // Set for response time calculation
        insertEvent(list, departureEvent);
        sim->totalServiceTime += serviceTime;
    } else {
        // Queue the customer
        if (sim->maxQueueSize == -1 || (sim->currentQueueLength + sim->inService) < sim->maxQueueSize) {
            sim->currentQueueLength++;
            if (sim->waitingCount >= sim->waitingCapacity) {
                sim->waitingCapacity = sim->waitingCapacity > 0 ? sim->waitingCapacity * 2 : 1;
                sim->waitingTimes = realloc(sim->waitingTimes, sim->waitingCapacity * sizeof(double));
            }
            sim->waitingTimes[sim->waitingCount++] = event->time; // Store arrival time for waiting calculation
        } else {
            // Arrival is dropped if the queue is at capacity
            sim->customerDropped++;
        }
    }

    // Schedule next arrival if needed
    if (sim->customerServed + sim->customerDropped < sim->totalCustomers) {
        double nextArrivalTime = sim->currentTime + exponentialRandom(sim->arrivalRate);
        insertEvent(list, createEvent(nextArrivalTime, ARRIVAL));
    }
}

void processDeparture(Simulation* sim, Event* event, EventList* list) {
    // Update queue length accounting for time passed since the last event
    double timeElapsed = event->time - sim->lastEventTime;
    sim->totalQueueLength += sim->currentQueueLength * timeElapsed;

    sim->currentTime = event->time;
    sim->customerServed++;
    sim->inService--;

    // Handle waiting time for the departing customer
    double responseTime = event->time - event->arrivalTime;
    sim->totalResponseTime += responseTime;

    if (sim->currentQueueLength > 0) {
        // Calculate waiting time for the next customer to start service
        double waitingTime = sim->currentTime - sim->waitingTimes[0];
        sim->totalWaitingTime += waitingTime;

        // Shift waiting times
        for (int i = 1; i < sim->waitingCount; i++) {
            sim->waitingTimes[i - 1] = sim->waitingTimes[i];
        }
        sim->waitingCount--;

        // Serve the next customer
        sim->currentQueueLength--;
        sim->inService++;
        double serviceTime = exponentialRandom(sim->serviceRate);
        Event* nextDepartureEvent = createEvent(sim->currentTime + serviceTime, DEPARTURE);
        if (sim->waitingCount > 0) {
            nextDepartureEvent->arrivalTime = sim->waitingTimes[0]; // Use the waiting time of the next in line
        } else {
            nextDepartureEvent->arrivalTime = sim->currentTime; // For cases without a waiting queue
        }
        insertEvent(list, nextDepartureEvent);
        sim->totalServiceTime += serviceTime;
    } else {
        // If no one is waiting and no service is happening, the system becomes idle
        sim->idleStartTime = sim->currentTime; // Mark the start of an idle period
    }
}

void computedStatistics(Simulation* sim) {
    double lambda = sim->arrivalRate;
    double mu = sim->serviceRate;
    int m = sim->servers;
    double rho = lambda / (mu * m);
    double p0 = 0.0;
    double En = 0.0, Er = 0.0, Ew = 0.0;
    double lambda_eff;

    // Handle infinite capacity case
	if(sim->maxQueueSize == -1 && lambda >= m * mu){
        sim->computedEn = HUGE_VAL;
        sim->computedEr = HUGE_VAL;
        sim->computedEw = HUGE_VAL;
        sim->computedP0 = 0.0;
        sim->computedLambdaEff = lambda; // The effective arrival rate remains lambda
	} else if (sim->maxQueueSize == -1) {
        // M/M/m formulae apply
        double sum = 0.0;
        for (int n = 0; n < m; ++n) {
            sum += pow(m * rho, n) / tgamma(n + 1);
        }
        double extraTerm = pow(m * rho, m) / (tgamma(m + 1) * (1 - rho));
        sum += extraTerm;
        p0 = 1 / sum;

        double omega = (pow(m * rho, m) * p0) / (tgamma(m + 1) * (1 - rho));
        En = m * rho + (rho * omega) / (1 - rho);
        lambda_eff = lambda; // No blocking, so effective arrival rate is actual arrival rate
        Er = En / lambda_eff;
        Ew = omega / (mu * (1 - rho));
		sim->computedEn = En;
		sim->computedEr = Er;
		sim->computedEw = Ew;
		sim->computedP0 = p0;
		sim->computedLambdaEff = lambda_eff;
    } else {
        double sum = 0.0;
        for (int n = 0; n <= sim->maxQueueSize; ++n) {
            if (n < m) {
                sum += pow(m * rho, n) / tgamma(n + 1);
            } else {
                sum += pow(m * rho, m) * pow(rho, n - m) / tgamma(m + 1);
            }
        }
        p0 = 1 / sum;

        double pB = (pow(m * rho, m) * pow(rho, sim->maxQueueSize - m)) / (tgamma(m + 1) * sum);
        lambda_eff = lambda * (1 - pB);

        for (int n = 0; n <= sim->maxQueueSize; ++n) {
            if (n < m) {
                En += n * (pow(m * rho, n) / tgamma(n + 1)) / sum;
            } else {
                En += n * (pow(m * rho, m) * pow(rho, n - m) / tgamma(m + 1)) / sum;
            }
        }
        Ew = (En / lambda_eff) - (1 / mu);
        Er = Ew + (1 / mu);
		sim->computedEn = En;
		sim->computedEr = Er;
		sim->computedEw = Ew;
		sim->computedP0 = p0;
		sim->computedLambdaEff = lambda_eff;
    }

}

void printStatistics(Simulation* sim, int departureCount, int finalPrint) {
    if (!finalPrint) {
        printf("\nAfter %d departures\n", departureCount);
    } else {
        printf("\nEnd of Simulation - after %d departures\n", departureCount);
    }

    // Check if the computedEn is infinity or not a number
    if (isinf(sim->computedEn) || isnan(sim->computedEn)) {
        printf("Mean n = %.4f (Simulated) and inf (Computed)\n", sim->totalQueueLength / sim->currentTime);
    } else {
        printf("Mean n = %.4f (Simulated) and %.4f (Computed)\n", sim->totalQueueLength / sim->currentTime, sim->computedEn);
    }

    // Repeat checks for computedEr and computedEw
    if (isinf(sim->computedEr) || isnan(sim->computedEr)) {
        printf("Mean r = %.4f (Simulated) and inf (Computed)\n", sim->totalResponseTime / sim->customerServed);
    } else {
        printf("Mean r = %.4f (Simulated) and %.4f (Computed)\n", sim->totalResponseTime / sim->customerServed, sim->computedEr);
    }

    if (isinf(sim->computedEw) || isnan(sim->computedEw)) {
        printf("Mean w = %.4f (Simulated) and inf (Computed)\n", sim->totalWaitingTime / sim->customerServed);
    } else {
        printf("Mean w = %.4f (Simulated) and %.4f (Computed)\n", sim->totalWaitingTime / sim->customerServed, sim->computedEw);
    }

    // For p0 and Effective arrival rate, infinity check is not needed
    printf("p0 = %.4f (Simulated) and %.4f (Computed)\n", sim->idleTime / sim->currentTime, sim->computedP0);
    printf("Effective arrival rate = %.4f (Simulated) and %.4f (Computed)\n", (double)sim->customerServed / sim->currentTime, sim->computedLambdaEff);
}


int main(int argc, char* argv[]) {
    if (argc < 8) {
        printf("Insufficient number of arguments provided!\n");
        return 1;
    }

    double lambda = atof(argv[1]);
    double mu = atof(argv[2]);
    int m = atoi(argv[3]);
    int B = atoi(argv[4]);
    int P = atoi(argv[5]); // Period for statistics printing
    int D = atoi(argv[6]); // Total number of customers to serve or total departures
    int seed = atoi(argv[7]);

	// Check for negative or zero values in the input parameters
    if (m <= 0 || D <= 0 || P <= 0) {
        printf("Input Error. Terminating Simulation...\n");
        return 1;
    }
    printf("Simulating M/M/m/B queue with lambda = %.6f, mu = %.6f, m = %d, B = %d, P = %d, D = %d, S = %d\n", lambda, mu, m, B, P, D, seed);

    Simulation* sim = initializeSimulation(lambda, mu, m, B, D, seed);
    EventList* eventList = createEventList();

    // Initialize the first arrival event
    double firstArrivalTime = exponentialRandom(sim->arrivalRate);
    insertEvent(eventList, createEvent(firstArrivalTime, ARRIVAL));

    int departureCount = 0;

    while (departureCount < D) {
        if (eventList->head == NULL) {
            fprintf(stderr, "Event list is empty before reaching desired number of departures.\n");
            break;
        }

        Event* event = eventList->head;
        eventList->head = event->next; // Remove the event from the list

        // Update current time to the event's time
        sim->currentTime = event->time;

        if (event->type == ARRIVAL) {
            processArrival(sim, event, eventList);
        } else if (event->type == DEPARTURE) {
            processDeparture(sim, event, eventList);
            departureCount++;
            
            // Print statistics periodically
            if (departureCount % P == 0 && departureCount != D) {
                computedStatistics(sim); // Ensure computed stats are updated before printing
                printStatistics(sim, departureCount, 0);
            }
        }

        free(event); // Free the processed event
    }

	printf("\ntotal_departures=%d  departure_count=%d  dropped_arrivals=%d\n", D, departureCount, sim->customerDropped);

    // Final print including computed statistics
    computedStatistics(sim); // Update computed statistics before final print
    printStatistics(sim, departureCount, 1); // Final print indicating simulation end

    // Free allocated resources
    while (eventList->head != NULL) {
        Event* temp = eventList->head;
        eventList->head = eventList->head->next;
        free(temp);
    }
    free(eventList);
	free(sim->waitingTimes);
    free(sim);
	

    return 0;
}
