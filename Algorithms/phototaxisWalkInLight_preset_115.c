// Written by Matan
// Updated 20210309
// Main changes
// 1. Stall communication for 10 minutes allowing reward systemt to reach steady state.

// Main prevoius changes:
// 1. Introduce finite tumble rate in light TUMBLERATEINLIGHT (helps avoid walls)
// 2. Flash state (green/blue) instead of keeping it on.

// INPUTS:
// kilo_turn_right - when 0, no communication
// kilo_straight_left - presets w0. (when 255, choose randomly from list)
// kilo_straight_right - presets w1. (when 255, choose randomly from list)

//All light intensity measurements are normalized to 0-255 (to be persistent with message)
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include "kilolib.h"

//#define DEBUG
//#include "debug.h"

//CHANGE IT BACK to 1:

#define TIMECONST 1 //Sets the time constant of the whole system. 1 is real time. Use 0.1 for rapid debugging.
#define TPS 32*TIMECONST //Ticks per second (see kilo_ticks)


//Light measurement params

#define MAXLIGHT 1023 //maximum value of 10bit light sensor 

#define QUICKLIGHTSAMPLES 100 //how many times to sample the light in a single measurement (used to be 100)

#define LIGHTSAMPLES 2 //How many persistent times to measure light in the integrator (stabilizes readings)

#define LIGHTSAMPLEINTERVAL 0.25*TPS

#define RELMARGINALREWARD 75 //minimal difference between rewards for exchange

//Communication Params
#define DISTANCETHRESHOLD 250 // How close robots need to be in order to exchange a message (55 is about touching)
//MOTION  PARAMS

//Motion-action states.
#define STAND   0
#define RUN     1
#define TUMBLE  2
#define ORBIT   3

//Motor values
#define STANDMOTVAL 0 //default motor val for standing (usually 0)
#define RUNMOTVAL 130  //default motor val for running Fast
#define ORBITMOTVALLOW 90 //70  //default motor val for running Slow (when light is ON)
#define ORBITMOTVALHIGH 130 //150  //default motor val for running Slow (when light is ON)
#define TUMBLEMOTVAL 110  //default motor val for tumbling (apply only to one motor)

//Motor action DURATIONS
#define STANDDURATION   30*TPS
#define RUNDURATION     2*TPS
#define TUMBLEDURATION  0.75*TPS
#define ORBITDURATION   6*TPS

#define NTUMBLEDURATIONS 6

#define TUMBLERATEINLIGHT 0.1//how often does a robot tumble when it thinks it's in the light (1/val).

#define COMMDELAY 600*TPS // Equivalent to about 10 minutes, roughly the time it takes to reach steady state.

#define MUTATIONRATE 0 //value between 0 (no mutation) to 1 (mutate every time) by currentDelay, and randomly picks a new strategy from the pool of strategies

// *****Neuron Temperature  ****
#define TNEURON 0.1
//Value between 0 and 1 that sets how much noise the action potential has. 
//When 0, no noise and neuron shoots only when >0.5
//when 1, full noise and neuron shoots when > rand_soft()/255;
//T is the error rate for shooting (sometimes, helps being persistently wrong).


// **** Globals Initionalization  *****
//Measured Light variables

double lightIntegral=0;
uint8_t currentLight = 255/2; //measured of current light intensity.
uint8_t lightSamples=0, lightFlag=1; //Light Flag is used to mark change in light status to initiate change in motor vals


//preset of different light state Thresholds. Dont forget to update NLIGHTSTATES to the length of this array.
#define NLIGHTSTATES 4

 
const uint8_t w0Presets[NLIGHTSTATES] = {1,4,16,64};//,100,120,140,160}; 
//(lower w1 means higher threshold)
const uint8_t w1Presets[NLIGHTSTATES] = {2,8,32,128};//,108,130,150,170}; //{30,150}; //number of different light threshold kilobots can have.
//The magic numbers above were chosen as pairs that exactly switch at the light region.

char runTumbleFlag = RUN; //takes values STAND, RUN, TUMBLE, ORBIT
int currentDelay = 0.1*TPS;

//Perceptron Params
//The perceptron is defined as a logistic function: P = 1/(1+e^(-(x0 w0+x1 w1))
//And normalized so that  give rise to values between close to 0 and close to 1 for P
//given that x1, w0,w1 are between 0 and 255.

double bias = -255;//this will be constant throughout (will be used for neuron)
//uint8_t w0=75, w1=128; //The weights for the perceptron.
uint8_t w0=75, w1=128; //The weights for the perceptron.

#define LIGHTMEM 250// Buffer size for light intensity measurements (light memory)
uint8_t lightMem[LIGHTMEM] = {0}; //same type as message; start with 0, initialize randomly later.
uint8_t lightMem1[LIGHTMEM] = {0}; //same type as message; start with 0, initialize randomly later.

uint8_t reward=0; //This is the reward function. Typical implementation is averaging the lightMem array.

//messaging variables

uint8_t distance_threshold=55;
message_t msg_tx;
// Flag to keep track of message transmission.
int message_sent = 0;

// Flag to keep track of new messages.
int new_message = 0;
uint8_t commOn=1; //this is used to turn on communication between bots. Set kilo_turn_right= 0 to supress comm.

//this mimics a long tailed distribution that for some reason the compiler struggles with:
//don't forget to update NTUMBLEDURATIONS when chaning array's length.
uint16_t tumbleDurations[] = {1*TUMBLEDURATION, 1*TUMBLEDURATION, 1.5*TUMBLEDURATION, 5*TUMBLEDURATION, 6*TUMBLEDURATION,7*TUMBLEDURATION};

//update the message to be transmitted
void updateMessageTx()
{
    if(commOn)
    {
        msg_tx.data[0] = reward;
        msg_tx.data[1] = w0;
        msg_tx.data[2] = w1;
        msg_tx.crc = message_crc(&msg_tx);
    }
    
    // It's important that the CRC is computed after the data has been set;
    // otherwise it would be wrong and the message would be dropped by the
    // receiver.
}

void message_rx(message_t *msg_rx, distance_measurement_t *distance)
{
    if (commOn) //commOn flags wheather communication is used. Set kilo_turn_right=0 to suppress messaging
    {
        // Set the flag on message reception.
        //uint8_t dist = estimate_distance(distance);
        if (kilo_ticks > COMMDELAY){ //Stall message processing until enough time has passed (~10minutes)
			uint8_t dist = estimate_distance(distance);
			if (dist<DISTANCETHRESHOLD){
				new_message = 1; //may be depricated
				//set_color(RGB(3,3,3));
				//delay(50);
				//set_color(RGB(0,0,0));
				
				//accept new threshold if sender's reward function is higher by a margin (10)
				
				if( (reward + RELMARGINALREWARD < msg_rx->data[0]))// && ((lightIntegral/lightSamples/1.)/RELMARGINALREWARD <= msg_rx->data[0]) )
				{
					w0 =  msg_rx->data[1];
					w1 =  msg_rx->data[2];
					updateMessageTx();
					set_color(RGB(1,3,3));
					delay(15);
					set_color(RGB(0,0,0)); // Light bulb effect (I learned!)
					delay(20);
					set_color(RGB(1,3,3));
					delay(15);
					set_color(RGB(0,0,0)); // Light bulb effect (I learned!)
				}
				
			}   
		}
    }
}

message_t *message_tx()
{
    return &msg_tx;
}

void message_tx_success()
{
    // Set the flag on message transmission.
    message_sent = 1;
}

void initMessaging()
{
    if (commOn) //use kilo_turn_right=0 to suppress messaging
    {
        // Register the message_tx callback function.
        kilo_message_tx = message_tx;
        // Register the message_rx callback function.
        kilo_message_rx = message_rx;
        // Register the message_tx_success callback function.
        kilo_message_tx_success = message_tx_success;
    }
}


double perceptron()//uint8_t w0, uint8_t w1)
{
    double x0,x1;
    double h,p;
	
    x0 = bias;
    x1 = lightIntegral*(1./lightSamples/1.);
    
//uint16_t bias = -255;//this will be constant throughout (will be used for neuron)
//uint8_t w0=75, w1=128; //The weights for the perceptron.
    //w0=75;
    x0=-255.;
    //w1=128;
    
    h = (w0*x0+w1*x1); 
    //h = (((double)w0)*x0+((double)w1)*x1)/1023.; //the 1023. prefactor is to normlize response within uint_8 range
    //h = (-255*75+128*x1)/1023.; 
    p = 1./(1.+exp(-1.*h));

	//if (rand_soft()/255.<TNEURON) //Introduce errors at a rate T
	//	p = 1-p;
    return p; //note that bias is always negative (haha) bias=-255.
    //return tanh(w0*bias+w1*currentLight);
}

uint8_t normalizeLightIntensity(uint32_t light)
{
    return (255*light)/MAXLIGHT;
}
//Measure light intensity at current position. Update global currentLight.
void sampleLight()
{
    // The ambient light sensor gives noisy readings. To mitigate this,
    // we take the average of 30 samples in quick succession. (be careful here. Too many samples and the int16 overflows)
    
    int numberOfSamples = 0;
    uint32_t sum = 0;
    int sample;
    while (numberOfSamples < QUICKLIGHTSAMPLES)
    {
        sample = get_ambientlight();
        
        // -1 indicates a failed sample, which should be discarded.
        if (sample != -1)
        {
            sample = get_ambientlight(); //normalize light intensity
            sum = sum + sample;
            numberOfSamples++;
        }
    }

    // Compute the average.
    currentLight = normalizeLightIntensity(sum / numberOfSamples);
}

//Sample light multiple times to get a spatio-temporal stable reading. Update global lightIntegral.
void integrateLight()
{
    //static unsigned long int lastTick=0;
    static uint32_t lastTick=0;

    if(kilo_ticks>lastTick + LIGHTSAMPLEINTERVAL) //light sampling intervals 
    {
        sampleLight();
        
		//The following little conditional allows lightIntegral never to vanish.
		if(lightSamples >= LIGHTSAMPLES){
			lightSamples =0;					
			lightIntegral = currentLight;
		}else
			lightIntegral += currentLight;
		
		lightSamples++;
                
        lastTick = kilo_ticks;
    }
}

//Insert integrated light measrument into light memory in a modular buffer (rolling buffer)
void updateLightBuffer()
{
    static int bufferState = 0; 
    if(bufferState == 0){
        set_color(RGB(0,3,0)); //Flash green when buffer re-sets  
        delay(50);
        set_color(RGB(0,0,0));
    } 
    bufferState = bufferState%LIGHTMEM; //make buffer state modular (keep rolling)
   
    lightMem[bufferState] = lightIntegral/lightSamples; // Insert integrated light into current circular buffer position
    
    bufferState++;  //advance buffer state.
}


void computeReward()
{
    int i;
    uint32_t sum=0;
    
    sum = 0;
    for(i=0;i<LIGHTMEM;i++)
        sum+=lightMem[i];
    
    reward = sum/LIGHTMEM;

    //fprintf(ptf, "REWARD : %d\n", reward);
}


double neuralThreshold()
{
    return (1.-TNEURON)*0.5 + (rand_soft()/255.)*TNEURON;
}
//Make an action based on light measuremnet (LED status and Persitent time update)
void indicateLight()
{
    if(lightSamples>=LIGHTSAMPLES)
    {
        
        updateLightBuffer(); //Store the measured light intensity
        computeReward();
        updateMessageTx(); //update the transmitted message.
        
        
        if( perceptron() >= 0.5)//neuralThreshold()) //much light, very bright.
        {
            if ((runTumbleFlag != RUN) & (runTumbleFlag!=STAND)){ //If not ALREADY doing one of these, then change to do so.
                runTumbleFlag = (rand_soft()%2)?RUN:STAND;    //here a random decision between RUN and STAND might be good
				currentDelay = 0; //Get out of current state.
				//runAndTumble();
			}
            set_color(RGB(0,1,0));
            delay(50);
            set_color(RGB(0,0,0));
        }
        else 
        {
            if ((runTumbleFlag != TUMBLE) & (runTumbleFlag!=ORBIT)){ //If not ALREADY doing one of these, then change to do so.
                runTumbleFlag = (rand_soft()%2)?TUMBLE:ORBIT;
				currentDelay = 0; //Get out of current state.
				//runAndTumble();
			}
            
           set_color(RGB(0,0,3));
           delay(50);
           set_color(RGB(0,0,0));
        }
        
        lightFlag = 1;
        //lightIntegral = 0; //This was a bad move as it sets to 0 every round.
        //lightSamples = 0;
    }
}


void initArray(uint8_t *arr, uint8_t arrSize, uint8_t val)
{
    static uint8_t i;
    
    for(i=0;i<arrSize;i++)
    {
        arr[i] = val;
    }
    
}
void setupMessaging()
{
    if(commOn)
    {
        msg_tx.type = NORMAL;
        // Some dummy data as an example.
        updateMessageTx();
    
        // Override message sending and cancel if kilo_turn_right is set to zero
        msg_tx.type = (commOn)?msg_tx.type:(-1);  
    }
}

void runAndTumble()
{
    //static unsigned long int lastTick=0;
    static uint32_t lastTick = 0;
    static uint8_t mleft,mright; //these will be used to set the motors values.
    
    //static uint8_t explore = 0; // this is a counter and every 100 steps the robots switches to an explorative tumble/orbit mode for some time.
    
    if((kilo_ticks > lastTick + currentDelay))// || lightFlag) //this means it's time to change direction!
    {
        switch ( runTumbleFlag )
        {
            case RUN:
                currentDelay = RUNDURATION; //Set the next delay to be that of the
                runTumbleFlag = STAND;  //Switch to stand for next change in mot val
                mleft = RUNMOTVAL;
                mright = RUNMOTVAL;
				
				
                break;
            case STAND:
                currentDelay  = STANDDURATION; //Set the next delay to be that of the
                runTumbleFlag = RUN;    //Switch to run for next change in mot val
                mleft = STANDMOTVAL;
                mright = STANDMOTVAL;
                
                break;
            
            case TUMBLE:
                currentDelay  = tumbleDurations[rand_soft()%NTUMBLEDURATIONS]; //Set the next delay to be that of the
                runTumbleFlag = ORBIT;
                
                if(rand_soft()%2) //turn left
                { 
                    mleft = TUMBLEMOTVAL;
                    mright = 0;
                }    
                else
                {
                    mleft = 0;
                    mright = TUMBLEMOTVAL;
                }
                break;
            
            case ORBIT:
                currentDelay = ORBITDURATION; //Set the next delay to be that of the
                runTumbleFlag = TUMBLE;
                
                
                switch (rand_soft()%3)
                {
                    case 0: //orbit CCW
                        mleft = ORBITMOTVALHIGH;
                        mright = ORBITMOTVALLOW;
                        break;
                    case 1: //orbit CW
                        mleft = ORBITMOTVALLOW;
                        mright = ORBITMOTVALHIGH;
                        break;
                    case 2:
                        mleft = ORBITMOTVALHIGH;
                        mright = ORBITMOTVALHIGH;
                        break;
                }
            default: ;
        }

        spinup_motors();
        set_motors(mleft,mright);
        //set_color(RGB(0,1,0));
        lightFlag = 0;
        lastTick = kilo_ticks;
        
    }
}

void mutate()
{
    //static unsigned long int lastTick=0;
    static uint32_t lastTick=0;

    if(kilo_ticks>lastTick + currentDelay) //light sampling intervals
    {
        if(rand_soft()<255*MUTATIONRATE){
            w0 = w0Presets[rand_soft()%NLIGHTSTATES];
            w1 = w1Presets[rand_soft()%NLIGHTSTATES];
        }
            
        lastTick = kilo_ticks;
    }
}

void setup() {
    //PLEASE add random delay to desynchronize
    
    //debug_init();
    
    commOn = 0; //use kilo_straight_right = 0 to supress messaging
    initMessaging();
    
    //randomly choose a light state from the different predefined thresholds (or set using kilo_straight_xxx)
    
    w0 = 1;
    w1 = 2;
    
    
    lightSamples = LIGHTSAMPLES;
    //Make an initial light measurment.
    sampleLight();
    
    //Extrapulate as though same intensity has been measured for some time:
    lightIntegral = currentLight*lightSamples;
    //initialize light memory buffer to the first measured value
    initArray(lightMem, LIGHTMEM, currentLight);
    updateLightBuffer();
    
    //Compute initial reward:
    computeReward();
    
    //Update message with reward and light threshold:
    updateMessageTx();
    
    set_color(RGB(1, 0, 0));
    delay(9500+rand_soft()*4); //delay for 20 seconds +-0.5
}

void loop() {
    runAndTumble();
    integrateLight();
    //mutate();
    indicateLight();
}

int main() {
    // initialize hardware
    kilo_init();
   
    // start program
    kilo_start(setup, loop);

    return 0;
}


