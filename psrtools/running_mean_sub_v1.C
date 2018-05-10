/***************************************************************************
 *
 *   Copyright (C) 2018 by Vishal Gajjar
 *   Licensed under the Academic Free License version 2.1
 *   A modified version of Paul's running_mean_sub to remove running mean from archive files. 
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
};


running_mean_sub::running_mean_sub ()
  : Application ("running_mean_sub", "running mean subtraction")
{
    do_weight = false;
}


void running_mean_sub::add_options (CommandLine::Menu& menu)
{
    CommandLine::Argument* arg;
    arg = menu.add (do_weight, 'w');
    arg->set_help("apply scaling to weights rather than to data");
}

void running_mean_sub::setup ()
{
}


void running_mean_sub::process (Pulsar::Archive* archive)
{
    // Window size
    const unsigned nbin = archive->get_nbin();		

    int nblock = 5; // Number of blocks across profile for smoothing 
    int win = floor(nbin/nblock);
    int extra = nbin - nblock*win;	
    cout << " Window " << win << " Blocks " << nblock << endl;
    cout << " extra " << extra << endl;
    	
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
		    //int iblock = 0;
		    //for (unsigned ibin=0; ibin<nbin;ibin++){    
		     for (int iblock = 0; iblock<nblock; iblock++){      
			    	float meandata = 0;
				for(int w=0;w<win | iblock+win>nbin ;w++)
					meandata += data[iblock*win+w];			   
				meandata /= win;
				for(int w=0;w<win;w++){
					data[iblock*win+w] /= meandata;							
					//printf("%d ",iblock*win+w);
				}
				// For the extra bins, normalize by the mean of the previous block
				if(nbin - ((iblock+1)*win) == extra)
					for(int e=0;e<extra;e++){
						data[(iblock+1)*win+e] /= meandata;
						//printf("%d ",(iblock+1)*win+e);
					}
	    	     }			
		    //}
		    // Write out a full chan
		    archive->get_Profile (isub, ipol, ichan)->set_amps(data);
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

