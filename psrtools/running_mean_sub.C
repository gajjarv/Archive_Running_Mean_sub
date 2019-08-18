/***************************************************************************
 *
 *   Copyright (C) 2018 by Vishal Gajjar
 *   Licensed under the Academic Free License version 2.1
 *   A modified version of Paul's normlize_rms to remove running mean from archive files. 
 *   Only useful for extracted files (non folded profiles, FRB type)
 ***************************************************************************/

using namespace std;

#include "Pulsar/Application.h"
#include "Pulsar/StandardOptions.h"
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/ProfileShiftFit.h"
#include "Pulsar/PhaseWeight.h"
#include "strutil.h"
#include "stdio.h"

#include "gsl/gsl_multifit.h"
#include "polifitgsl.h"

using namespace Pulsar;

//
//! Archive running mean subtractor 
//
class running_mean_sub : public Pulsar::Application
{
public:

    //! Default constructor
    running_mean_sub ();

    //! Setup
    void setup ();

    //! Process the given archive
    void process (Pulsar::Archive*);

protected:
    void add_options (CommandLine::Menu&);

    bool do_weight;
    int deg;
};


running_mean_sub::running_mean_sub ()
  : Application ("running_mean_sub", "running mean subtraction")
{
    do_weight = false;
    deg = 2;
}


void running_mean_sub::add_options (CommandLine::Menu& menu)
{
    CommandLine::Argument* arg;
    //arg = menu.add (do_weight, 'w');
    //arg->set_help("apply scaling to weights rather than to data");
    arg = menu.add(deg,"d");
    // This is not working yet
    arg->set_help("Order of polynomial fit (not working yet, auto calc now)");
}

void running_mean_sub::setup ()
{
}

void running_mean_sub::process (Pulsar::Archive* archive)
{
    
    int nbin = archive->get_nbin();		
    if(nbin>4000) deg = 10;
    if(nbin<4000 and nbin>2000) deg = 8;
    if(nbin<2000 and nbin>1000) deg = 6;
    if(nbin<1000 and nbin>150) deg = 5;    
    if(nbin<150) deg = 2;
    int nblock = 20; // Number of blocks across profile for smoothing 
    int win = floor(nbin/nblock);
    int extra = nbin - nblock*win;	
    // cout << " Window " << win << " Blocks " << nblock << endl;
    // cout << " extra " << extra << endl;
    	
    // Output filename
    cout << "# " << archive->get_filename() << endl;

    Reference::To<Archive> copy = archive->clone();
    copy->tscrunch();
    copy->convert_state(Signal::Intensity);

    bool scale_indep_subints = true;

    // Get baseline values
    // Use this for applying single scaling vs freq to all subints
    vector< vector< Estimate<double> > > mean;
    vector< vector< double> > var;
    if (scale_indep_subints == false)
        copy->get_Integration(0)->baseline_stats(&mean, &var);
	
    printf("Polynomial fit with degree : %d\n",deg);	
    // Dedisperse makes no difference as we are operating on individule channels
   // archive->dedisperse(); 	
    // Loop over subints
    for (unsigned isub=0; isub<archive->get_nsubint(); isub++) {

        // Use this to get a new scaling for each subint
        if (scale_indep_subints)
            archive->get_Integration(isub)->baseline_stats(&mean, &var);

        // Loop over channels
        for (unsigned ichan=0; ichan<archive->get_nchan(); ichan++) {
            double scale;
            if (var[0][ichan]>0.0) {
                //scale = 1.0 / sqrt(var[0][ichan]); // Single-pol or stokes
		scale=1.0;
                for (unsigned ipol=0; ipol<archive->get_npol(); ipol++) 
                {
		    //printf("Scale: %lf \n",scale);
		    //scale=1.0;
                    //archive->get_Profile(isub,ipol,ichan)->scale(scale);
		
		    // Read a full chan 
		    float* data = archive->get_Profile (isub, ipol, ichan)->get_amps();
		    
		    // FIRST Technique
		    /* Using blocks. The extracted time-series gets divided in to many blocks 
		     * and average are determined from each block for subtraction (didn't work very well)
		     *
		     int iblock = 0;
		     for (int iblock = 0; iblock<nblock; iblock++){      
			    	float meandata = 0;
				for(int w=0;w<win | iblock+win>nbin ;w++)
					meandata += data[iblock*win+w];			   
				meandata /= win;
				for(int w=0;w<win;w++){
					data[iblock*win+w] -= meandata;							
					//printf("%d ",iblock*win+w);
				}
				// For the extra bins, normalize by the mean of the previous block
				if(nbin - ((iblock+1)*win) == extra)
					for(int e=0;e<extra;e++)
						data[(iblock+1)*win+e] -= meandata;
	    	     }
		     */		

		    // SECOND Technique
		    /* Moving aveage is calculated for each bin and subtracted. 
		     * (worked but too many free parameters to adjust)
		     
		    float oldmeandata=0;
		    for (unsigned ibin=0; ibin<nbin;ibin++){
			float meandata=0;   
			if(ibin+win>=nbin){
				for(; ibin<nbin;ibin++){
					data[ibin] -= oldmeandata;			
				}
			}	
			else{
				for(int w=0; w<win ; w++){
                         	       meandata += data[ibin+w];
                        	}
				meandata /= win;
				data[ibin] -= meandata; 
			}	
			oldmeandata=meandata;
		    }
		    archive->get_Profile (isub, ipol, ichan)->set_amps(data);
		    */
		    //

		    // THIRD Technique
		    
		    float* x; 
		    float* smooth;
		    //int deg=30;
		    double coeff[deg];
		    x = (float *)malloc(nbin * sizeof(float));
		    smooth = (float *)malloc(nbin * sizeof(float));
		    
		    for (unsigned ibin=0; ibin<nbin;ibin++) x[ibin] = ibin;
		    
	    	    smooth = polynomialfit(nbin, deg, x, data, coeff);	    

		    // Write out a full chan
		    archive->get_Profile (isub, ipol, ichan)->set_amps(smooth);
		    
		    //
                }
            } else {
                archive->get_Integration(isub)->set_weight(ichan, 0.0);
            }
            // Reset all non-zero weights to 1
            // XXX why did I do this.... ?
#if 0 
            if (archive->get_Integration(isub)->get_weight(ichan)!=0.0) 
                archive->get_Integration(isub)->set_weight(ichan,1.0);
#endif
          
        }
    }

    // Unload corrected archive
    // TODO use standard unload options
    archive->unload(replace_extension(archive->get_filename(), "norm"));

}

static running_mean_sub program;

int main (int argc, char** argv)
{
    return program.main (argc, argv);
}

