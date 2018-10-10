#ifndef HistoTools_cxx
#define HistoTools_cxx

// User
#include "../interface/HistoTools.h"

// ROOT
#include <TFitResult.h>
#include <TF1.h>
#include <TH1.h>

//#define DEBUG

//****************************************************************************
void HistoTools::BookHisto_1D(TH1D *&histo, 
			      const char *hName,
			      const char *hTitle, 
			      const Int_t nBins, 
			      const Double_t xMin,  
			      const Double_t xMax)
//****************************************************************************
{

  
  histo = new TH1D(hName, hTitle, nBins, xMin, xMax);
  histo->Sumw2();
  return;
}


//****************************************************************************
void HistoTools::BookHisto_2D(TH2D *&histo, 
			      const char *hName,
			      const char *hTitle, 
			      const Int_t nBinsX, 
			      const Double_t xMin,  
			      const Double_t xMax,
			      const Int_t nBinsY, 
			      const Double_t yMin,  
			      const Double_t yMax)
//****************************************************************************
{

  histo = new TH2D(hName, hTitle, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  histo->Sumw2();
  return;
}


//****************************************************************************
void HistoTools::DivideHistos_1D(TH1D *hNumerator,
				 TH1D *hDenominator)
//****************************************************************************
{

   // Sanity Check
   for (int i = 0; i <= hNumerator->GetNbinsX()+1; i++){
     int N    = hNumerator  ->GetBinContent(i);
     int D    = hDenominator->GetBinContent(i);
     double r = double(N)/double(D);
     if ( r > 1.0) std::cout << "=== HistoTools::DivideHistos_1D()  - bin " << i 
   			    << " (eT = " << i*hNumerator->GetBinWidth(i) << " GeV): Numerator/Denominator = " 
   			    << N << "/" << D << " = " << r << std::endl;
   }
   
   hNumerator->Divide(hNumerator, hDenominator, 1.0, 1.0, "binomial");
  return;
}


//****************************************************************************
TFitResultPtr HistoTools::FitFunctionSL(TH1D *h,
					TF1 *fitFunction,
					Option_t* fitOptions,
					Option_t* graphicOptions,
					std::vector<double> v_fitArea,
					double significanceLevel,
					bool bPrintFitInfo,
					bool bSaveFitInfo)
 //****************************************************************************
 {
		     
   // Variable declaration
   Double_t fitRangeLow  = -1.0;
   Double_t fitRangeHigh = -1.0;
   Double_t fitArea   = -1.0;
   Double_t chiSq     = -1.0;
   Double_t dof       = -1.0;
   Double_t redChiSq  = +99999.9;
   Double_t area      = -1.0;
   Double_t areaError = -1.0;
   Double_t mean      = -1.0;
   Double_t meanError = -1.0;
   Double_t sigma     = -1.0;
   Double_t sigmaError= -1.0;
   bool bFitSuccess   = false;
   //
   Double_t chiSq_tmp          = -1.0;
   Double_t chiSq_quantileLow  = -1.0;
   Double_t chiSq_quantileHigh = -1.0;
   Double_t dof_tmp            = -1.0;
   Double_t redChiSq_tmp       = -1.0;
   Double_t area_tmp           = -1.0;
   Double_t areaError_tmp      = -1.0;
   Double_t mean_tmp           = -1.0; 
   Double_t meanError_tmp      = -1.0;
   Double_t sigma_tmp          = -1.0; 
   Double_t sigmaError_tmp     = -1.0;
   bool bFitSuccess_tmp        = false;
   // Probability of rejecting hyporthesis. The smaller the alpha the wider the allowed RedChiSq range.
   // 95% Confidence Level is alpha=0.025: chiSq_{alpha}=0.975 to chiSq_{alpha}=0.025
   // Christiana (22 Oct 2015): To significance level (type I error sou, oti theleis na to apokaleseis) en to alpha in total.
   // ara eheis miso alpha pou ti mia k miso alpha pou tin alli ara eheis sta deksia alpha/2 kai sta aristera 1-(alpha/2)
   double alpha      = significanceLevel/2.0;
   double alpha_High = 1.0 - (alpha);
   double alpha_Low  = alpha;

   // Create a table with all the fit information
   Table fitInfo("Histogram | Area | xFit-Low | xFit-High | Chi^2 | DOF | RedChi^2 | Norm. | Mean | Sigma | RedChi^2 Range | Fit ", "Text");

   int counter = 0;
   // Loop over all fit-area values and attempt to fit in the corresponding symmetric fit range
   for(std::vector<double>::iterator it = v_fitArea.begin(); it != v_fitArea.end(); ++it, counter++) {
  
     // Get the area to look for fit range
     Double_t fitArea_tmp = *it;

     // Find symmetric range to perform the fit
     Double_t fitRangeLow_tmp  = 0.0;
     Double_t fitRangeHigh_tmp = 0.0;
     FindSymmetricFitRange(h, fitArea_tmp, fitRangeLow_tmp, fitRangeHigh_tmp);

     // Get fit results to given fit range but do not plot the fit result ("0")
     TString fitOpts(fitOptions);
     fitOpts = fitOpts + "0";

     TFitResultPtr r_tmp = FitFunctionToRange(h, fitFunction, fitOpts.Data(), graphicOptions, fitRangeLow_tmp, fitRangeHigh_tmp);
     if (!r_tmp->IsEmpty() && !r_tmp->IsZombie() ){
     chiSq_tmp          = r_tmp->Chi2();
     dof_tmp            = r_tmp->Ndf();
     redChiSq_tmp       = r_tmp->Chi2()/dof_tmp;
     area_tmp           = r_tmp->Parameter(0);
     areaError_tmp      = r_tmp->ParError(0);
     mean_tmp           = r_tmp->Parameter(1); 
     meanError_tmp      = r_tmp->ParError(1); 
     sigma_tmp          = r_tmp->Parameter(2); 
     sigmaError_tmp     = r_tmp->ParError(2);
     chiSq_quantileLow  = TMath::ChisquareQuantile(alpha_High, dof_tmp);  // for lower limit on redChiSq
     chiSq_quantileHigh = TMath::ChisquareQuantile(alpha_Low , dof_tmp);  // for upper limit on redChiSq
     }
  
     // Check that fit meets certain minimum criteria. Rule of thumb:
     // 0) a chiSq = 1 indicates the perfect fit.
     // 1) a chiSq >> 1 indicates a poor model fit,
     // 2) a chiSq >  1 indicates that the fit has not fully captured the data (or that the error variance has been underestimated).
     // 3) A chiSq <  1 indicates that the model is 'over-fitting' the data: either the model is improperly fitting noise, or the error variance has been overestimated
     // Evaluate the quantiles of the chi-squared pdf for a given probability value (p) and degrees of freedom (dof)
     // Quantile alpha is the area (probability) to the right of the chi^2 value. The quantile 1-alpha is the probability
     // to the left of the chi^2 value (i.e. up to and including that value)
     chiSq_quantileLow  = TMath::ChisquareQuantile(alpha_Low , dof_tmp) / dof_tmp;
     chiSq_quantileHigh = TMath::ChisquareQuantile(alpha_High, dof_tmp) / dof_tmp;
     if( (redChiSq_tmp >= chiSq_quantileLow) && (redChiSq_tmp <= chiSq_quantileHigh) ) bFitSuccess_tmp = true;
     else bFitSuccess_tmp = false;
     
     // Update table with fit properties
     fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitArea_tmp) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitRangeLow_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitRangeHigh_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( chiSq_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( dof_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( redChiSq_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( area )  + " +/- " + tools.ToString( areaError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( mean )  + " +/- " + tools.ToString( meanError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( sigma ) + " +/- " + tools.ToString( sigmaError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( chiSq_quantileLow ) + " - " + tools.ToString( chiSq_quantileHigh ) );
     fitInfo.AddRowColumn(counter, tools.ToString( bFitSuccess_tmp ) );

     // Check that fit is a success and even better than before
     if (!bFitSuccess_tmp) continue;
     if ( fabs(1.0-redChiSq_tmp) >= fabs(1.0-redChiSq) ) continue;
    
     // Update fit results only if new redChiSq is better than last.
     fitArea      = fitArea_tmp;
     fitRangeLow  = fitRangeLow_tmp;
     fitRangeHigh = fitRangeHigh_tmp;
     chiSq        = chiSq_tmp;
     dof          = dof_tmp;
     redChiSq     = redChiSq_tmp;
     area         = area_tmp;
     areaError    = areaError_tmp;
     mean         = mean_tmp;
     meanError    = meanError_tmp;
     sigma        = sigma_tmp;
     sigmaError   = sigmaError_tmp;
     bFitSuccess  = bFitSuccess_tmp;
  
   }//eof: for(std::vector<double>::iterator it = v_fitArea.begin(); it != v_fitArea.end(); ++it, counter++) {

   // Re-draw the best fit only if fit is success
   TString fitOptsFinal(fitOptions);
   if(!bFitSuccess) fitOptsFinal = fitOptsFinal + "N0";
   else fitOptsFinal = fitOptsFinal + "";
   TFitResultPtr r = FitFunctionToRange(h, fitFunction, fitOptsFinal.Data(), graphicOptions, fitRangeLow, fitRangeHigh);

   // Extend the fit range (indicate that it is not part of the fit by line style and width)
   TF1 *f =  (TF1*) fitFunction->Clone("fitFunction_clone");
   f->SetLineWidth(2);
   f->SetLineStyle(kDashed);
   f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
   if(bFitSuccess) h->GetListOfFunctions()->Add(f); 

   // Save Fit info to file or print?
   fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitArea) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitRangeLow ) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitRangeHigh ) );
   fitInfo.AddRowColumn(counter, tools.ToString( chiSq ) );
   fitInfo.AddRowColumn(counter, tools.ToString( dof ) );
   fitInfo.AddRowColumn(counter, tools.ToString( redChiSq ) );
   fitInfo.AddRowColumn(counter, tools.ToString( area )  + " +/- " + tools.ToString( areaError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( mean )  + " +/- " + tools.ToString( meanError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( sigma ) + " +/- " + tools.ToString( sigmaError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( chiSq_quantileLow ) + " - " + tools.ToString( chiSq_quantileHigh ) );
   fitInfo.AddRowColumn(counter, tools.ToString( bFitSuccess ) );
   if(bPrintFitInfo) fitInfo.Print();
   if(bSaveFitInfo) fitInfo.SaveToFile("fitInfo.txt", "a");
   
   return r;
 }



//****************************************************************************
TFitResultPtr HistoTools::FitFunction(TH1D *h,
				      TF1* fitFunction,
				      Option_t* fitOptions,
				      Option_t* graphicOptions,
				      double fitRangeLow,
				      double fitRangeHigh,
				      std::vector<double> v_fitArea,
				      double redChiSqMax,
				      double redChiSqMin)
//****************************************************************************
 {
		     
   // Optionally fit over the entire histogram range
   if (fitRangeLow == -1.0 && fitRangeHigh == -1.0){
     fitRangeLow  = h->GetXaxis()->GetBinLowEdge(0);
     fitRangeHigh = h->GetXaxis()->GetBinLowEdge( h->GetNbinsX()+1 );
   }

   // Variable declaration
   Double_t fitArea   = -1.0;
   Double_t chiSq     = -1.0;
   Double_t dof       = -1.0;
   Double_t redChiSq  = +99999.9;
   Double_t area      = -1.0;
   Double_t areaError = -1.0;
   Double_t mean      = -1.0;
   Double_t meanError = -1.0;
   Double_t sigma     = -1.0;
   Double_t sigmaError= -1.0;
   bool bFitSuccess   = false;
   //
   Double_t chiSq_tmp          = -1.0;
   Double_t chiSq_quantileLow  = -1.0;
   Double_t chiSq_quantileHigh = -1.0;
   Double_t dof_tmp            = -1.0;
   Double_t redChiSq_tmp       = -1.0;
   Double_t area_tmp           = -1.0;
   Double_t areaError_tmp      = -1.0;
   Double_t mean_tmp           = -1.0; 
   Double_t meanError_tmp      = -1.0;
   Double_t sigma_tmp          = -1.0; 
   Double_t sigmaError_tmp     = -1.0;
   bool bFitSuccess_tmp        = false;

   // Create a table with all the fit information
   Table fitInfo("Histo Name | Fit Area | xFit-Low | xFit-High | Chi^2 | DOF | RedChi^2 | Norm. | Norm. Error | Mean | Mean Error | Sigma | Sigma Error | Fit | Chosen", "Text");
 
   // Loop over all fit-area values and attemptto fit in the corresponding symmetric fit range
   int counter = 0;
   for(std::vector<double>::iterator it = v_fitArea.begin(); it != v_fitArea.end(); ++it, counter++) {
  
     // Get the area to look for fit range
     Double_t fitArea_tmp = *it;

     // Find symmetric range to perform the fit
     Double_t fitRangeLow_tmp  = 0.0;
     Double_t fitRangeHigh_tmp = 0.0;
     FindSymmetricFitRange(h, fitArea_tmp, fitRangeLow_tmp, fitRangeHigh_tmp);

     // Get fit results to given fit range but do not plot the fit result ("0")
     TString fitOpts(fitOptions);
     fitOpts = fitOpts + "0";

     TFitResultPtr r_tmp = FitFunctionToRange(h, fitFunction, fitOpts.Data(), graphicOptions, fitRangeLow_tmp, fitRangeHigh_tmp);
     if (!r_tmp->IsEmpty() && !r_tmp->IsZombie() ){
     chiSq_tmp          = r_tmp->Chi2();
     dof_tmp            = r_tmp->Ndf();
     redChiSq_tmp       = r_tmp->Chi2()/dof_tmp;
     area_tmp           = r_tmp->Parameter(0);
     areaError_tmp      = r_tmp->ParError(0);
     mean_tmp           = r_tmp->Parameter(1); 
     meanError_tmp      = r_tmp->ParError(1); 
     sigma_tmp          = r_tmp->Parameter(2); 
     sigmaError_tmp     = r_tmp->ParError(2);
     }
  
     // Check that fit meets certain minimum criteria. Rule of thumb:
     // 0) a chiSq = 1 indicates the perfect fit.
     // 1) a chiSq >> 1 indicates a poor model fit,
     // 2) a chiSq >  1 indicates that the fit has not fully captured the data (or that the error variance has been underestimated).
     // 3) A chiSq <  1 indicates that the model is 'over-fitting' the data: either the model is improperly fitting noise, or the error variance has been overestimated
     if( (dof_tmp > 2) && (dof_tmp < 150)  && (redChiSq_tmp <= redChiSqMax) && (redChiSq_tmp >= redChiSqMin) ) bFitSuccess_tmp = true;
     else bFitSuccess_tmp = false;
      
     // Update table with fit properties 
     fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitArea_tmp) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitRangeLow_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( fitRangeHigh_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( chiSq_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( dof_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( redChiSq_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( area ) );
     fitInfo.AddRowColumn(counter, tools.ToString( areaError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( mean ) );
     fitInfo.AddRowColumn(counter, tools.ToString( meanError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( sigma ) );
     fitInfo.AddRowColumn(counter, tools.ToString( sigmaError ) );
     fitInfo.AddRowColumn(counter, tools.ToString( chiSq_quantileLow ) + " - " + tools.ToString( chiSq_quantileHigh ) );
     fitInfo.AddRowColumn(counter, tools.ToString( bFitSuccess_tmp ) );
     fitInfo.AddRowColumn(counter, tools.ToString( false ) );

     // Check that fit is a success and even better than before
     if (!bFitSuccess_tmp) continue;
     if ( fabs(1.0-redChiSq_tmp) >= fabs(1.0-redChiSq) ) continue;
    
     // Update fit results only if new redChiSq is better than last.
     fitArea      = fitArea_tmp;
     fitRangeLow  = fitRangeLow_tmp;
     fitRangeHigh = fitRangeHigh_tmp;
     chiSq        = chiSq_tmp;
     dof          = dof_tmp;
     redChiSq     = redChiSq_tmp;
     area         = area_tmp;
     areaError    = areaError_tmp;
     mean         = mean_tmp;
     meanError    = meanError_tmp;
     sigma        = sigma_tmp;
     sigmaError   = sigmaError_tmp;
     bFitSuccess  = bFitSuccess_tmp;
  
   }//eof: for(std::vector<double>::iterator it = v_fitArea.begin(); it != v_fitArea.end(); ++it, counter++) {

   // Re-draw the best fit only if fit is success
   TString fitOptsFinal(fitOptions);
   if(!bFitSuccess) fitOptsFinal = fitOptsFinal + "N0";
   else fitOptsFinal = fitOptsFinal + "";
   TFitResultPtr r = FitFunctionToRange(h, fitFunction, fitOptsFinal.Data(), graphicOptions, fitRangeLow, fitRangeHigh);

   // Extend the fit range (indicate that it is not part of the fit by line style and width)
   TF1 *f =  (TF1*) fitFunction->Clone("fitFunction_clone");
   f->SetLineWidth(2);
   f->SetLineStyle(kDashed);
   f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
   if(bFitSuccess) h->GetListOfFunctions()->Add(f); 

   // Save Fit info to file or print?
   fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitArea) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitRangeLow ) );
   fitInfo.AddRowColumn(counter, tools.ToString( fitRangeHigh ) );
   fitInfo.AddRowColumn(counter, tools.ToString( chiSq ) );
   fitInfo.AddRowColumn(counter, tools.ToString( dof ) );
   fitInfo.AddRowColumn(counter, tools.ToString( redChiSq ) );
   fitInfo.AddRowColumn(counter, tools.ToString( area ) );
   fitInfo.AddRowColumn(counter, tools.ToString( areaError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( mean ) );
   fitInfo.AddRowColumn(counter, tools.ToString( meanError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( sigma ) );
   fitInfo.AddRowColumn(counter, tools.ToString( sigmaError ) );
   fitInfo.AddRowColumn(counter, tools.ToString( bFitSuccess ) );
   fitInfo.AddRowColumn(counter, tools.ToString( true ) );
   fitInfo.Print();
   fitInfo.SaveToFile("fitInfo.txt", "a");
   
   return r;
 }


//****************************************************************************
TFitResultPtr HistoTools::FitFunctionAlt(TH1D *h,
				      TF1 *fitFunction,
				      Option_t* fitOptions,
				      Option_t* graphicOptions,
				      double fitRangeLow,
				      double fitRangeHigh,
				      std::vector<double> v_fitArea,
				      double redChiSqMax,
				      double redChiSqMin)
//****************************************************************************
{
  
  // 5-sigma iterative outlier removal guassian fit alternative. Preliminary stages. So-and-so results
  
  // Create a table with all the fit information
  Table fitInfo("Histo Name | Fit Area | xFit-Low | xFit-High | Chi^2 | DOF | RedChi^2 | Norm.  | Mean | Sigma | Fit Success | Chosen", "Text");

  // Variable declaration
  v_fitArea.clear();
  std::vector<Double_t> v_chiSq;
  std::vector<Double_t> v_dof;
  std::vector<Double_t> v_redChiSq;
  std::vector<Double_t> v_area;
  std::vector<Double_t> v_areaError;
  std::vector<Double_t> v_mean;
  std::vector<Double_t> v_meanError;
  std::vector<Double_t> v_sigma;
  std::vector<Double_t> v_sigmaError;
  std::vector<Double_t> v_fitRangeLow;
  std::vector<Double_t> v_fitRangeHigh;
  std::vector<bool> v_bFitSuccess;

  //
  TString fitOpts(fitOptions);
  fitOpts = fitOpts + "0";
  fitRangeLow  = h->GetXaxis()->GetBinLowEdge(0);
  fitRangeHigh = h->GetXaxis()->GetBinLowEdge( h->GetNbinsX()+1 );
  TFitResultPtr r = FitFunctionToRange(h, fitFunction, fitOpts.Data(), graphicOptions, fitRangeLow, fitRangeHigh);
  Double_t fitArea   = 1.0;
  Double_t chiSq     = r->Chi2();
  Double_t dof       = r->Ndf();
  Double_t redChiSq  = r->Chi2()/dof;
  Double_t area      = r->Parameter(0);
  Double_t areaError = r->ParError(0);
  Double_t mean      = r->Parameter(1); 
  Double_t meanError = r->ParError(1); 
  Double_t sigma     = r->Parameter(2); 
  Double_t sigmaError= r->ParError(2); 

  // Start with initial hit on entire histogram range
  bool bFitSuccess = ( (redChiSq <= redChiSqMax) && (redChiSq >= redChiSqMin) );
  if(bFitSuccess){
    v_fitArea.push_back(fitArea);
    v_chiSq.push_back(chiSq);
    v_dof.push_back(dof);
    v_redChiSq.push_back(redChiSq);
    v_area.push_back(area);
    v_areaError.push_back(areaError );
    v_mean.push_back(mean);
    v_meanError.push_back(meanError);
    v_sigma.push_back(sigma);
    v_sigmaError.push_back(sigmaError);
    v_fitRangeLow.push_back(fitRangeLow);
    v_fitRangeHigh.push_back(fitRangeHigh);
    v_bFitSuccess.push_back(bFitSuccess);
  }

  
  // For-loop: Over all histogram bins, removing outliers
  int counter = 0;
  for( int i=1; i < h->GetNbinsX()+1-1; i++){

    // Check that for this fit range there are no 5-sigma outliers
    if( ( fabs(fitRangeLow) < (mean + 5*sigma) ) || fabs(fitRangeHigh) < (mean + 5*sigma) ) continue;

    // Calculate new fit range
    fitRangeLow  = h->GetXaxis()->GetBinCenter(i);
    fitRangeHigh = h->GetXaxis()->GetBinCenter( h->GetNbinsX()+1-i );
    fitArea      = (h->Integral( i, h->GetNbinsX()+1-i )) / h->Integral( 0, h->GetNbinsX()+1 );

    // Ensure that fit encloses at least 50% of histogram
    if (fitArea < 0.5) break;

    // Re-fit in a range that doesn't include the 5-sigma outliers
    r = FitFunctionToRange(h, fitFunction, fitOpts.Data(), graphicOptions, fitRangeLow, fitRangeHigh);
    chiSq     = r->Chi2();
    dof       = r->Ndf();
    redChiSq  = r->Chi2()/dof;
    area      = r->Parameter(0);
    areaError = r->ParError(0);
    mean      = r->Parameter(1); 
    meanError = r->ParError(1); 
    sigma     = r->Parameter(2); 
    sigmaError= r->ParError(2);
      
    // Check that fit meets certain minimum criteria
    bFitSuccess = ( (redChiSq <= redChiSqMax) && (redChiSq >= redChiSqMin) );
    if(!bFitSuccess) continue;
    
    // Updated fit values
    v_fitArea.push_back(fitArea);
    v_chiSq.push_back(chiSq);
    v_dof.push_back(dof);
    v_redChiSq.push_back(redChiSq);
    v_area.push_back(area);
    v_areaError.push_back(areaError );
    v_mean.push_back(mean);
    v_meanError.push_back(meanError);
    v_sigma.push_back(sigma);
    v_sigmaError.push_back(sigmaError);
    v_fitRangeLow.push_back(fitRangeLow);
    v_fitRangeHigh.push_back(fitRangeHigh);
    v_bFitSuccess.push_back(bFitSuccess);

    // Update table with fit properties
    fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
    fitInfo.AddRowColumn(counter, tools.ToString( fitArea) );
    fitInfo.AddRowColumn(counter, tools.ToString( fitRangeLow ) );
    fitInfo.AddRowColumn(counter, tools.ToString( fitRangeHigh ) );
    fitInfo.AddRowColumn(counter, tools.ToString( chiSq ) );
    fitInfo.AddRowColumn(counter, tools.ToString( dof ) );
    fitInfo.AddRowColumn(counter, tools.ToString( redChiSq ) );
    fitInfo.AddRowColumn(counter, tools.ToString( area )  + " +/- " + tools.ToString( areaError ) );
    fitInfo.AddRowColumn(counter, tools.ToString( mean )  + " +/- " + tools.ToString( meanError ) );
    fitInfo.AddRowColumn(counter, tools.ToString( sigma ) + " +/- " + tools.ToString( sigmaError ) );
    fitInfo.AddRowColumn(counter, tools.ToString( bFitSuccess ) );
    fitInfo.AddRowColumn(counter, tools.ToString( false ) );
    counter++;

  }//eof:for( int i=0; i < h->GetNbinsX()+1; i++){
  

  // Any successful fits?
  std::vector<bool>::iterator it = std::find(v_bFitSuccess.begin(), v_bFitSuccess.end(), true);
  if (it == v_bFitSuccess.end() ){
    FitFunctionToRange(h, fitFunction, fitOpts.Data(), graphicOptions, 1, -1);
    return r;
  }
 
  
  int iBestFit = 0;
  // For-loop: All successful fits
  for(int j=0; j < (int) v_redChiSq.size(); j++){

    // Only consider successful fits
    if(!v_bFitSuccess.at(j)) continue;

    // Find reduced chi^2 closest to 1
    double OneMinusRedChiSq_tmp = fabs( 1.0 - v_redChiSq.at(j) );
    double OneMinusRedChiSq     = fabs( 1.0 - v_redChiSq.at(iBestFit) );
    
    if ( OneMinusRedChiSq_tmp >= OneMinusRedChiSq ) continue;   
    else iBestFit=j;

  }//eof: for(int j=0; j < (int) v_redChiSq.size(); j++){

  iBestFit = v_dof.size()-1;
  
  // Re-draw the best fit only if fit is success
  TString fitOptsFinal(fitOptions);
  FitFunctionToRange(h, fitFunction, fitOptsFinal.Data(), graphicOptions, v_fitRangeLow.at(iBestFit), v_fitRangeHigh.at(iBestFit));
  
  // Extend the fit range (indicate that it is not part of the fit by line style and width)
  TF1 *f =  (TF1*) fitFunction->Clone("fitFunction_clone");
  f->SetLineWidth(2);
  f->SetLineStyle(kDashed);
  f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  h->GetListOfFunctions()->Add(f); 

  // Last table row is the actual fit 
  fitInfo.AddRowColumn(counter, tools.ToString( h->GetName() ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_fitArea.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_fitRangeLow.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_fitRangeHigh.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_chiSq.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_dof.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_redChiSq.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_area.at(iBestFit) )  + " +/- " + tools.ToString( v_areaError.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_mean.at(iBestFit) )  + " +/- " + tools.ToString( v_meanError .at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_sigma.at(iBestFit) ) + " +/- " + tools.ToString( v_sigmaError.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( v_bFitSuccess.at(iBestFit) ) );
  fitInfo.AddRowColumn(counter, tools.ToString( true ) );    
  fitInfo.Print();
  fitInfo.SaveToFile("fitInfo.txt", "a");

  
  return r;
}
 //****************************************************************************



 //****************************************************************************
 TFitResultPtr HistoTools::FitFunctionToRange(TH1D *h,
					      TF1 *fitFunction,
					      Option_t* fitOptions,
					      Option_t* graphicOptions,
					      double fitRangeLow,
					      double fitRangeHigh)
 //****************************************************************************
 {
		     
   // Fit the custom function to the histogram passed as argument   
   TFitResultPtr r = h->Fit(fitFunction, fitOptions, graphicOptions, fitRangeLow, fitRangeHigh);

   return r;
}


//****************************************************************************
void HistoTools::FindSymmetricFitRange(TH1D *h,
				       double fitArea,
				       double &fitRangeLow,
  				       double &fitRangeHigh,
				       bool bMakeSymmetric)
//****************************************************************************
{
  
  // Get histograms characteristics
  int binPeak       = h->GetMaximumBin();
  double xPeak      = h->GetXaxis()->GetBinCenter(binPeak);
  double normFactor = h->Integral( 0, h->GetNbinsX()+1 );

  // Loop from peak bin to the right until the integral equals the fitArea
  for( int i=binPeak; i < h->GetNbinsX()+1; i++){
    
    double binValue   = h->GetXaxis()->GetBinCenter(i);
    double binContent = h->GetBinContent(i);
    double coverage   = h->Integral(binPeak, i)/normFactor;
    if( coverage <= fitArea) fitRangeHigh  = h->GetXaxis()->GetBinCenter(i);
  }
  
  // Loop from peak bin to the left until the integral is 0.477 of normFactor
  for( int j=binPeak; j >= 0; j--){
    
    double binValue   = h->GetXaxis()->GetBinCenter(j);
    double binContent = h->GetBinContent(j);
    double coverage   = h->Integral(j, binPeak)/normFactor;
    if( coverage <= fitArea) fitRangeLow = h->GetXaxis()->GetBinCenter(j);
  }
  
  // Ensure symmetric fit range (Gives much better fits)
  if(!bMakeSymmetric) return;
  if ( fabs(fitRangeHigh) > fabs(fitRangeLow) ) fitRangeLow = -fitRangeHigh;
  else fitRangeHigh = -fitRangeLow;

  return;
}



//****************************************************************************
void HistoTools::ConvertToOneMinusCumulativeHisto_1D(TH1D *histo)
//****************************************************************************
{
  
  // For-loop over all histogram bins
  const int nBins = histo->GetNbinsX();
  
  for (int b=1; b <= nBins; b++){   //starting from bin 1 (because bin 0 is underflow)
    const double value = histo->Integral(b, nBins);
    const double error = TMath::Sqrt(value);
    histo->SetBinContent( b, value );
    histo->SetBinError  ( b, error );
  }
  
  return;
}


//****************************************************************************
void HistoTools::ConvertToOneMinusCumulativeHisto_2D(TH2D *histo)
//****************************************************************************
{
  // Leadng jet on x-axis

  // For-loop over all histogram bins
  const int nBinsX = histo->GetNbinsX();
  const int nBinsY = histo->GetNbinsY();
  
  // For-loop: x-axis bins
  for (int bx=1; bx <= nBinsX; bx++){

    // For-loop: y-axis bins
    for (int by=1; by <= nBinsY; by++){

      if (by > bx) break;
      const double value = histo->Integral(bx, nBinsX, by, nBinsY);
      const double error = TMath::Sqrt(value);
      histo->SetBinContent( bx, by, value );
      histo->SetBinError( bx, by, error );

    } // For-loop: y-axis bins
    
  } // For-loop: x-axis bins

  return;
}


//****************************************************************************
void HistoTools::ConvertToRateHisto_1D(TH1D *histo,
				       const Double_t nEvents, 
				       const Double_t crossingRate)
//****************************************************************************
{
  /*
    To normalize trigger rates, you can use the
    fact that your background sample (NeutrinoGun aka Minimum Bias) 
    corresponds to the crossing rate, since, as soon as PU >> 1, 
    every crossing leads to an interaction. The crossing rate is not 
    40 MHz because not all the bunches are filled. We usually take 
    30 MHz for our normalizations. 

    Hence, if you have n events that pass your trigger selection, out of
    the N events of your sample (e.g. N = 500k for the full NeutrinoGun
    sample at PU=140), then the rate is simply:
    Rate = 30 MHz x  (n/N)

    Note: default value for "crossingRate" is 30.0E+06
  */


  // Variable declarations
  const Double_t convertTokHz  = ( 1.0E-03);  // 1Hz -> 1kHz
  const Double_t normFactor    = (crossingRate * convertTokHz) / (nEvents);
  
    // Convert histogram to 1-cumulative
  ConvertToOneMinusCumulativeHisto_1D(histo);

  // Scale to desired crossing-rate
  histo->Scale( normFactor );


#ifdef DEBUG
  // std::cout << "nEvents = " << nEvents << std::endl;
  // std::cout << "histo->GetEntries() = " << histo->GetEntries() << std::endl;
  std::cout << "normFactor = Crossing_Rate x ConvertToKHz = " << crossingRate << " x " << convertTokHz << " = " << normFactor << std::endl;
#endif
  return;

}


//****************************************************************************
void HistoTools::ConvertToRateHisto_2D(TH2D *histo, 
				       const Double_t nEvents,
				       const Double_t crossingRate)
//****************************************************************************
{

  // Variable declarations
  const Double_t convertTokHz = ( 1.0E-03);  // 1Hz -> 1kHz
  const Double_t normFactor   = (crossingRate * convertTokHz) / (nEvents);

  // Convert histogram to 1-cumulative
  ConvertToOneMinusCumulativeHisto_2D(histo);

  // Scale to desired crossing-rate
  histo->Scale( normFactor );

  return;
}

//****************************************************************************
void HistoTools::FillAllBinsUpToValue_1D(TH1D *histo, 
					 const Double_t fillValue)
//****************************************************************************
{

  // Variable declarations
  const int stepValue = int ( double(histo->GetXaxis()->GetXmax())/double(histo->GetNbinsX()) ); 
  const int maxValue  = histo->GetXaxis()->GetXmax();

  // For-loop: All bins
  for (int value = 0; value <= maxValue; value += stepValue){ 
    if ( fillValue >= value) histo->Fill(value);
  }
  
  return;
}


//****************************************************************************
void HistoTools::FillAllBinsUpToValue_2D(TH2D *histo, 
					 const Double_t xFillValue, 
					 const Double_t yFillValue)
//****************************************************************************
{

  // Variable declarations
  const int stepValue = int ( double(histo->GetXaxis()->GetXmax())/double(histo->GetNbinsX()) ); 
  const int yMaxValue = histo->GetYaxis()->GetXmax();
  const int xMaxValue = histo->GetXaxis()->GetXmax();

  if (yMaxValue != xMaxValue)
    {
      std::cout << "=== HistoTools::FillAllBinsUpToCaloEtHisto_2D() - X-axis max is different than Y-axis max.\nExiting \n";
      exit(1);    
    }
  
  // For-loop: x-axis 
  for (int xVal = 0; xVal <= xMaxValue; xVal+=stepValue)
    {
      
      // Only x-axis values up to xFillValue
      if (xVal > xFillValue) break;
      
      // For-loop: y-axis
      for (int yVal = 0; yVal <= yMaxValue; yVal+=stepValue)
	{
	  
	  // Only y-axis values up to yFillValue
	  if (yVal > yFillValue) break;
	  
	  // Sanity check: xEt is leading value => xVal > yVal
	  if (yVal > xVal) break;
      
	  // Fill histogram with the values
	  histo->Fill(xVal, yVal);

	}// For-loop: y-axis 
      
    }// For-loop: x-axis
  
  return;
}


//****************************************************************************
void HistoTools::ConvertNumeratorHistoToEfficiency_1D(TH1D *histo, 
						      const int nEvtsTotal, 
						      const std::string errType)
//****************************************************************************
{
  
  // Variable declarations
  const int nBins = histo->GetNbinsX()+1;
  double eff, err;

  // For-loop: all bins
  for (int bx = 0; bx <= nBins; bx++){
    
    // Calculate efficiency value and error
    const int nPass = histo->GetBinContent(bx);
    tools.Efficiency(nPass, nEvtsTotal, errType, eff, err );
  
    // Replace histogram's bin contect with efficiency and its error
    histo->SetBinContent(bx, eff);
    histo->SetBinError  (bx, err);
  }

  return;
}


//****************************************************************************
void HistoTools::ConvertNumeratorHistoToEfficiency_2D(TH2D *histo, 
						      const int nEvtsTotal,
						      const std::string errType)
//****************************************************************************
{

  // Variable declarations
  const int nBinsX  = histo->GetNbinsX()+1;
  const int nBinsY  = histo->GetNbinsY()+1;
  double eff, err;
  
  // For-loop: x-axis bins
  for (int bx=0; bx <= nBinsX; bx++){

    // For-loop: y-axis bins
    for (int by=0; by <= nBinsY; by++){

      const int nPass = histo->GetBinContent(bx, by);
      tools.Efficiency(nPass, nEvtsTotal, errType, eff, err );

      // Update current histo bin to true eff value and error
      histo->SetBinContent(bx, by, eff);
      histo->SetBinError  (bx, by, err);

    } // For-loop: y-axis bins

  }// For-loop: x-axis bins

  return;
}

#endif  // HistoTools_cxx
