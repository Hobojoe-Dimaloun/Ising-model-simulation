/*********************
* Date of creation: 02/11/17
* Author: Michael O'Donnell
* Contact: mo14776@my.bristol.ac.uk
* Other Authors:
**************************************
* Change History
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

//Screen dimension constants
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 800;

//Starts up SDL and creates window
bool init();

//Loads media
bool loadMedia();

//Frees media and shuts down SDL
static void close();

void screenRenderfunc(int *v, int x, int y);

static void randLattice(int *lattice, int x, int y, int z);

//The window we'll be rendering to
SDL_Window* gWindow = NULL;

//The window renderer
SDL_Renderer* gRenderer = NULL;

const double gBoltzmann = 1.38064852E-23;
const double gPi = 3.14159;

int main(int argc,char **argv)
{
    int numOfNodes;
    int debug = 0;  // Debug
    int x = 800; // x - dimension lattice sites
    int y = 800; // y - dimension lattice sites
    int z = 1; // z - dimension latties sites
    double temperature = 0.0001; // Kelvin
    int J = 1;
    //
    // Handle input arguments
    //
    int numOfThreads;
    for(int i = 0; i < argc; i++)
    {
        switch (i)
        {
            case 0: numOfThreads = omp_get_max_threads(); break;
            case 1: numOfThreads =  atoi(argv[i]); break;
            case 2: numOfNodes = *argv[i]; break;
            case 3: debug = 1; break;
            default: break;

        }
    }

    omp_set_num_threads(numOfThreads);
    if(debug == 1)
    {
        printf("Warning! Debug mode entered. Press any key to continue...\n");
        getchar();
    }

    printf("Ising model begin.\n");

    int *ising_lattice = NULL;

    //
    // Allocate lattice memory
    //

    if(( ising_lattice = calloc(x * y * z, sizeof(*ising_lattice) )) == NULL)
    {
        printf("-Error %d : %s\nFile : %s\nLine : %d\n", errno, strerror(errno), __FILE__, __LINE__);
    }

    //
    // initialise SDL
    //
    if( !init() )
    {
        printf( "Failed to initialize!\n" );
        return 0;
    }

    //
    // Assign random spins to the lattice, ie 1 or -1
    //

    int loop = 0;

    randLattice(ising_lattice,x,y,z);

    gsl_rng *rndarray[numOfThreads];
    printf("%d\n", numOfThreads);
    for(int i = 0; i<numOfThreads; i++)
    {
        rndarray[i] = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(rndarray[i],i);
    }

    //
    //
    //
    bool quit = false;

	//Event handler
	SDL_Event e;
	//While application is running
	while( !quit )
	{


        #pragma omp parallel
        {
            gsl_rng *r = rndarray[omp_get_thread_num()];
            #pragma omp single
            {
    		//Handle events on queue
        		while( SDL_PollEvent( &e ) != 0 )
        		{
        			//User requests quit
        			if( e.type == SDL_QUIT)
        			{
        				quit = true;
        			}
                    if( e.type == SDL_KEYDOWN)
                    {
                        if( e.key.keysym.sym == SDLK_UP)
            			{
            				temperature+=0.01 ;
                            printf("T = %lf\n", temperature);

            			}
                        if( e.key.keysym.sym  == SDLK_DOWN)
            			{
            				temperature-=0.01 ;
                            printf("T = %lf\n", temperature);
            			}
                        if( e.key.keysym.sym  == SDLK_RIGHT)
            			{
                            randLattice(ising_lattice,x,y,z);
            			}
                        if( e.key.keysym.sym  == SDLK_LEFT)
            			{
                            J*=-1;
            			}

                    }

        		}


                if(loop%100 == 0)
                screenRenderfunc(ising_lattice,x,y);
                //
                // generate random point
                //
            }


            int location = (double)gsl_rng_uniform(r) * x * y ;
            int icolumn = location%y; // the remainder is the column number
            int irow = (location - icolumn)/x  ;// number of rows
            int temp;


            double dE = 0, E1=0, E2=0;
            //
            // Calculate energy of cell above
            //
            (irow == 0) ? (temp = x-1) : (temp = irow-1);

            E1+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell below
            //
            (irow == (x-1)) ? (temp = 0) : (temp = irow+1);

            E1+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell left
            //

            (icolumn == 0) ? (temp = y-1) : (temp = icolumn-1);
            E1+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Calculate energy of cell right
            //
            (icolumn == (y-1)) ? (temp = 0) : (temp = icolumn+1);
            E1+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Flip the spin
            //

            ising_lattice[location] *=-1 ;

            (irow == 0) ? (temp = x-1) : (temp = irow-1);

            E2+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell below
            //
            (irow == (x-1)) ? (temp = 0) : (temp = irow+1);

            E2+= J*ising_lattice[location] * ising_lattice[temp*x + icolumn]*(-1);

            //
            // Calculate energy of cell left
            //

            (icolumn == 0) ? (temp = y-1) : (temp = icolumn-1);
            E2+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            //
            // Calculate energy of cell right
            //
            (icolumn == (y-1)) ? (temp = 0) : (temp = icolumn+1);
            E2+= J*ising_lattice[location] * ising_lattice[irow*x + temp]*(-1);

            dE= E2 - E1;

            double U = gsl_rng_uniform(r);
            double V = gsl_rng_uniform(r);

            double z0 = sqrt(-2 * log(U)) * cos(2*gPi*V);
            double z1 = sqrt(-2 * log(U)) * sin(2*gPi*V);

            double guass = z1;
            //printf("guass %lf\n",guass );
            if( dE > 0)
            {

                double beta =1/temperature ;

                double prob = exp(beta *dE*(-1)) ;
            //    printf(" prob = %g\n",prob);

                if( prob < gsl_rng_uniform(r))
                {
                //    printf(" t = %g\n",dE);
                    ising_lattice[location]*= -1;

                }
            }

        }

        loop ++;
    }

    close();

    return 0;
}

void screenRenderfunc(int *v, int x, int y)
{
    const int pixelwidth = SCREEN_WIDTH /  x;
    const int pixelheight = SCREEN_HEIGHT /  y;

	//While application is running
	//Clear screen
	SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );
	SDL_RenderClear( gRenderer );

	//Render red filled quad

    for(int i = 0; i < x; i ++)
    {
        for (int j = 0; j < y; j ++)
        {
            SDL_Rect fillRect = { i*pixelwidth, j*pixelheight, pixelwidth,pixelheight};
			if(v[i*x +j] == 1)
            {
                SDL_SetRenderDrawColor( gRenderer, 0xFF, 0x00, 0x00, 0xFF );
            }
            else
            {
                SDL_SetRenderDrawColor( gRenderer, 0x00, 0xFF, 0x00, 0xFF );
            }
            SDL_RenderFillRect( gRenderer, &fillRect );

        }
    }

	//Update screen
	SDL_RenderPresent( gRenderer );
}

bool init()
{
	//Initialization flag
	bool success = true;

	//Initialize SDL
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 )
	{
		printf( "SDL could not initialize! SDL Error: %s\n", SDL_GetError() );
		success = false;
	}
	else
	{
		//Set texture filtering to linear
		if( !SDL_SetHint( SDL_HINT_RENDER_SCALE_QUALITY, "1" ) )
		{
			printf( "Warning: Linear texture filtering not enabled!" );
		}

		//Create window
		gWindow = SDL_CreateWindow( "SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
		if( gWindow == NULL )
		{
			printf( "Window could not be created! SDL Error: %s\n", SDL_GetError() );
			success = false;
		}
		else
		{
			//Create renderer for window
			gRenderer = SDL_CreateRenderer( gWindow, -1, SDL_RENDERER_SOFTWARE );
			if( gRenderer == NULL )
			{
				printf( "Renderer could not be created! SDL Error: %s\n", SDL_GetError() );
				success = false;
			}
			else
			{
				//Initialize renderer color
				SDL_SetRenderDrawColor( gRenderer, 0xFF, 0xFF, 0xFF, 0xFF );

				//Initialize PNG loading
				int imgFlags = IMG_INIT_PNG;
				if( !( IMG_Init( imgFlags ) & imgFlags ) )
				{
					printf( "SDL_image could not initialize! SDL_image Error: %s\n", IMG_GetError() );
					success = false;
				}
			}
		}
	}

	return success;
}

static void close()
{
	//Destroy window
	SDL_DestroyRenderer( gRenderer );
	SDL_DestroyWindow( gWindow );
	gWindow = NULL;
	gRenderer = NULL;

	//Quit SDL subsystems
	IMG_Quit();
	SDL_Quit();
}

static void randLattice(int *lattice, int x, int y, int z)
{
    for(int i = 0; i < x * y * z; i ++)
    {
        if((rand()/(double)RAND_MAX) <0.5)
        {
            lattice[i] = 1;
        }
        else
        {
            lattice[i] = -1;
        }
    }
}
