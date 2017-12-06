#include <iostream>
#include <vector>
#include <string>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "DataFormat/IOManager.h"
#include "DataFormat/EventImage2D.h"



int main( int nargs, char** argv ) {

  gStyle->SetOptStat(0);
  
  std::cout << "Check Instance image completeness" << std::endl;

  std::string inputlarcvfile = "larcv.root";
  if ( nargs==2 )
    inputlarcvfile = argv[1];
  
  larcv::IOManager io( larcv::IOManager::kREAD );
  io.add_in_file( inputlarcvfile );
  io.initialize();

  size_t nentries = io.get_n_entries();
  float maxadc = 80.0;
  bool dump_images = false;
  
  TFile* out = new TFile("output_check_instanceimg.root", "recreate" );
  
  // per-event tree
  TTree* evtree  = new TTree("evcheck", "Event-level label quality metrics" );
  int nnonlabeled[4]; // [U,V,Y,total]
  int nabovethresh[4];
  double fracmissed[4];

  int ancestor_nnonlabeled[4]; // [U,V,Y,total]
  int ancestor_nabovethresh[4];
  double ancestor_fracmissed[4];
  
  int nparticles;
  evtree->Branch( "nnonlabeled", nnonlabeled, "nnonlabeled[4]/I" );
  evtree->Branch( "nabovethresh", nabovethresh, "nabovethresh[4]/I" );
  evtree->Branch( "fracmissed", fracmissed, "fracmissed[4]/D" );

  evtree->Branch( "ancestor_nnonlabeled",  ancestor_nnonlabeled,  "ancestor_nnonlabeled[4]/I" );
  evtree->Branch( "ancestor_nabovethresh", ancestor_nabovethresh, "ancestor_nabovethresh[4]/I" );
  evtree->Branch( "ancestor_fracmissed",   ancestor_fracmissed,   "ancestor_fracmissed[4]/D" );
  
  evtree->Branch( "nparticles",  &nparticles, "nparticles/I" );
  
  
  for (size_t ientry=0; ientry<nentries; ientry++) {
    io.read_entry(ientry);
    std::cout << "====================" << std::endl;
    std::cout << "Entry " << ientry << std::endl;
    std::cout << "====================" << std::endl;

    
    larcv::EventImage2D* ev_instance = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, "instance" );
    larcv::EventImage2D* ev_ancestor = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, "ancestor" );
    larcv::EventImage2D* ev_wire     = (larcv::EventImage2D*)io.get_data( larcv::kProductImage2D, "wire" );
    
    const std::vector<larcv::Image2D>& instance_v = ev_instance->Image2DArray();
    const std::vector<larcv::Image2D>& ancestor_v = ev_ancestor->Image2DArray();    
    const std::vector<larcv::Image2D>& img_v      = ev_wire->Image2DArray();    

    std::cout << "Number of images in instance_v: " << instance_v.size() << std::endl;
    
    const larcv::ImageMeta& metap2 = instance_v.at(2).meta();
    std::set<int> trackidset;
    
    nnonlabeled[3]  = 0;
    nabovethresh[3] = 0;
    ancestor_nnonlabeled[3]  = 0;
    ancestor_nabovethresh[3] = 0;

    TCanvas clabel( "c", "", 3600, 1200 );
    clabel.Divide(3,2);

    TH2D* hplane[3]          = {NULL};
    TH2D* hplane_ancestor[3] = {NULL};    
    
    for ( size_t p=0; p<instance_v.size(); p++) {

      nnonlabeled[p]  = 0;
      nabovethresh[p] = 0;
      ancestor_nnonlabeled[p] = 0;
      ancestor_nabovethresh[p] = 0;
      
      const larcv::Image2D& labelimg    = instance_v[p];
      const larcv::Image2D& ancestorimg = ancestor_v[p];            
      const larcv::Image2D& img         = img_v[p];

      if ( dump_images ) {
	char pname[100];
	sprintf( pname, "hp%03d", (int)p );
	hplane[p]          = new TH2D( pname, "", metap2.cols(), metap2.min_x(), metap2.max_x(), metap2.rows(), metap2.min_y(), metap2.max_y() );
      
	char pname_ancestor[100];
	sprintf( pname_ancestor, "hp%03d_ancestor", (int)p );      
	hplane_ancestor[p] = new TH2D( pname_ancestor, "", metap2.cols(), metap2.min_x(), metap2.max_x(), metap2.rows(), metap2.min_y(), metap2.max_y() );
      }
      
      for (size_t irow=0; irow<metap2.rows(); irow++) {
	for (size_t icol=0; icol<metap2.cols(); icol++) {

	  float adc = img.pixel(irow,icol);
	  
	  if ( adc<10.0 )
	    continue; // below threshold

	  float filladc = adc;
	  if ( filladc>maxadc )
	    filladc = maxadc;

	  if ( dump_images ) {
	    hplane[p]->SetBinContent( icol+1, metap2.rows()-1-irow, filladc );
	    hplane_ancestor[p]->SetBinContent( icol+1, metap2.rows()-1-irow, filladc );
	  }

	  nabovethresh[p]++;
	  nabovethresh[3]++;
	  
	  int trackid    = (int)labelimg.pixel( irow, icol );
	  if ( trackid==-1.0 ) {
	    nnonlabeled[p]++;
	    nnonlabeled[3]++;
	    if ( dump_images )
	      hplane[p]->SetBinContent( icol+1, metap2.rows()-1-irow, -10 );
	  }
	  else {
	    trackidset.insert(trackid);
	  }
	  
	  int ancestorid = (int)ancestorimg.pixel( irow, icol );
	  if ( ancestorid==-1.0 ) {
	    ancestor_nnonlabeled[p]++;
	    ancestor_nnonlabeled[3]++;
	    if ( dump_images )
	      hplane_ancestor[p]->SetBinContent( icol+1, metap2.rows()-1-irow, -10 );
	  }
	  
	}
      }
      fracmissed[p] = double(nnonlabeled[p])/double(nabovethresh[p]);
      ancestor_fracmissed[p] = double(ancestor_nnonlabeled[p])/double(nabovethresh[p]);

      if ( dump_images ) {
	clabel.cd(p+1);      
	hplane[p]->SetMinimum(-11);
	hplane[p]->SetMaximum( maxadc);
	hplane[p]->Draw("COLZ");
	
	clabel.cd(3+p+1);      
	hplane_ancestor[p]->SetMinimum(-11);
	hplane_ancestor[p]->SetMaximum( maxadc);
	hplane_ancestor[p]->Draw("COLZ");
      }
      
    }//end of plane loop

    fracmissed[3]          = double(nnonlabeled[3])/double(nabovethresh[3]);
    ancestor_fracmissed[3] = double(ancestor_nnonlabeled[3])/double(nabovethresh[3]);    
    nparticles = trackidset.size();

    if ( dump_images ) {
      char imgname[100];
      sprintf( imgname, "imgid_ientry%03d.png", (int)ientry );
      clabel.Update();
      clabel.SaveAs( imgname );

      for (int p=0; p<3; p++) {
	delete hplane[p];
	delete hplane_ancestor[p];
      }
    }
    
    evtree->Fill();
    //break;
  }
  
  out->Write();
  
  return 0;
}
