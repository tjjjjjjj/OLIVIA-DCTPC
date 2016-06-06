#include "Deconvolution.hh" 


TH1 * waveform::deconv::WienerDeconv( const TH1 * signal, const waveform::FFT * H, const double * SNR, double noise_estimate) 
{

  //compute signal fft
  waveform::FFT sig_fft(signal->GetNbinsX()); 
  double * sig_vals = new double[signal->GetNbinsX()]; 
  fftw_complex * filter = new fftw_complex[signal->GetNbinsX()/2 + 1]; 

  TH1 * output = (TH1*) signal->Clone(TString(signal->GetName()) + TString("_deconvolved")); 

  
  for (int i = 0; i <signal->GetNbinsX(); i++)
  {
      sig_vals[i] = signal->GetBinContent(i+1); 
  }

  sig_fft.setWaveform(sig_vals); 
  sig_fft.fft(); 

  if (!SNR && noise_estimate < 0)
  {
    noise_estimate = sig_fft.getFFTMag2(signal->GetNbinsX()/2); 
  }

  double min = 1; 
  if (noise_estimate <=0)
  {
    for (int i = 0; i < signal->GetNbinsX(); i++)
    {
      if (sig_fft.getFFTMag2(i) < min && sig_fft.getFFTMag2(i) > 0)
      {
        min = sig_fft.getFFTMag2(i); 
      }

      noise_estimate = min*1e-6; //avoid singularity
    }
  }

  TH1F * filter_h = new TH1F("filter_h","filter_h",signal->GetNbinsX()/2 + 1, 0, signal->GetNbinsX()/2+1); 

  for (int i = 0; i <signal->GetNbinsX()/2 + 1; i++)
  {

      double snr = SNR ? SNR[i] : sig_fft.getFFTMag2(i) / noise_estimate; 
      if (snr < 0) snr = 0; 
//      std::cout << "snr: " << snr << std::endl; 
//      std::cout << "hreal: " << H->getFFTRealPart(i) << std::endl; 
//      std::cout << "hmag2: " << H->getFFTMag2(i) << std::endl; 
      filter[i][0]= H->getFFTRealPart(i) / (H->getFFTMag2(i) + 1./snr); 
      filter[i][1]= -H->getFFTImPart(i) / (H->getFFTMag2(i) + 1./snr); 
      std::cout << filter[i][0] << " + " << filter[i][1] << "i" << std::endl; 
      filter_h->SetBinContent(i+1, filter[i][0]*filter[i][0] + filter[i][1] * filter[i][1]);  
  }

  filter_h->Draw(); 

  sig_fft.multiply(filter); 
  sig_fft.inverse_fft();  

  for (int i = 0; i <signal->GetNbinsX(); i++)
  {
    std::cout << sig_fft.getInverseFFT(i) << std::endl; 
    output->SetBinContent(i+1,sig_fft.getInverseFFT(i)); 
  }

  return output; 

}


