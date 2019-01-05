#include <windows.h>
#include <fstream>
#include <iostream>
#include "mdk/mdk.h"
#include "math.h"

#include <complex>
#include <iostream>
#include <valarray>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef std::valarray<float> FArray;

#pragma optimize ("awy", on)

#define MACHINE_NAME "fFourierFx"
#define MACHINE_AUTHOR "Thor Muto Asmund"
#ifdef _DEBUG
#define MACHINE_VERSION "0.01.0 beta"
#define MACHINE_FULL_NAME "Nafai's fFourierFx (debug)"
#else
#define PROGRAM_VERSION "0.01.0"
#define MACHINE_FULL_NAME "Nafai's fFourierFx"
#endif

#define M_PI   3.1415926535897932384626433832795
#define M_2PI 6.283185307179586476925286766559005

#define WINDOWFN_SQUARE 0
#define WINDOWFN_HAMMING 1
#define WINDOWFN_HANN 2
#define WINDOWFN_GAUSS 3
#define WINDOWFN_FLAT_TOP 4
#define WINDOWFN_BLACKMAN 5
#define WINDOWFN_BLACKMAN_HARRIS 6
#define WINDOWFN_BLACKMAN_NUTTALL 7
#define WINDOWFN_NUTTALL 8
#define WINDOWFN_KAISER 9
#define WINDOWFN_KAISER_BESSEL 10
#define WINDOWFN_TRAPEZOID 11
#define WINDOWFN_SINC 12
#define WINDOWFN_SINE 13

#define OVERLAPCURVE_EXPONENTIAL 0
#define OVERLAPCURVE_LINEAR 1

void fft(CArray& x);
void ifft(CArray& x);
void fft_fast(CArray &x);
void ifft_fast(CArray& x);
float parameterValueToFloat(int param, word value);
int parameterValueToInt(int param, word value);
void create_window_data(FArray& window, int length, int windowType, double alpha, double beta, bool unityGain);

////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////

#define PARAM_FLOOR_IDX 0
#define PARAM_FLOOR_NONE 0xFF
#define PARAM_FLOOR_DEFAULT 0x00
CMachineParameter const paramFloor =
{
	pt_byte,                    // Parameter data type
	"Floor",                       // Parameter name as its shown in the parameter window
	"Noise floor",              // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	0xFE,                       // Maximum value
	PARAM_FLOOR_NONE,                // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_FLOOR_DEFAULT              // the default slider value
};

#define PARAM2_IDX 1
#define PARAM2_NONE 0xFFFF
#define PARAM2_DEFAULT 0x0000
CMachineParameter const param2 =
{
	pt_word,                    // Parameter data type
	"P2",                       // Parameter name as its shown in the parameter window
	"Parameter 2",              // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	0x8000,                     // Maximum value
	PARAM2_NONE,                // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM2_DEFAULT              // the default slider value
};

#define PARAM_BINSIZE_IDX 2
#define PARAM_BINSIZE_NONE 0xFF
#define PARAM_BINSIZE_DEFAULT 0x09 // 512
CMachineParameter const paramBinSize =
{
	pt_byte,                    // Parameter data type
	"Bin size",                 // Parameter name as its shown in the parameter window
	"FFT bin size",             // Parameter description as its shown in the pattern view's statusbar
	0x03,                       // Minimum value
	0x0F,                       // Maximum value
	PARAM_BINSIZE_NONE,         // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_BINSIZE_DEFAULT       // the default slider value
};

#define PARAM_WINDOWFN_IDX 3
#define PARAM_WINDOWFN_NONE 0xFF
#define PARAM_WINDOWFN_DEFAULT WINDOWFN_SQUARE
CMachineParameter const paramWindowingFunction =
{
	pt_byte,                    // Parameter data type
	"Window",                   // Parameter name as its shown in the parameter window
	"FFT windowing function",   // Parameter description as its shown in the pattern view's statusbar
	WINDOWFN_SQUARE,            // Minimum value
	WINDOWFN_SINE,              // Maximum value
	PARAM_WINDOWFN_NONE,        // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_WINDOWFN_DEFAULT      // the default slider value
};

#define PARAM_ALPHA_IDX 4
#define PARAM_ALPHA_NONE 0xFFFF
#define PARAM_ALPHA_DEFAULT 0
CMachineParameter const paramAlpha =
{
	pt_word,                    // Parameter data type
	"Alpha",                    // Parameter name as its shown in the parameter window
	"FFT windowing alpha",      // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	1000,                       // Maximum value
	PARAM_ALPHA_NONE,           // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_ALPHA_DEFAULT         // the default slider value
};

#define PARAM_BETA_IDX 5
#define PARAM_BETA_NONE 0xFFFF
#define PARAM_BETA_DEFAULT 0
CMachineParameter const paramBeta =
{
	pt_word,                    // Parameter data type
	"Beta",                     // Parameter name as its shown in the parameter window
	"FFT windowing beta",       // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	10000,                      // Maximum value
	PARAM_BETA_NONE,            // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_BETA_DEFAULT          // the default slider value
};

#define PARAM_OVERLAP_IDX 6
#define PARAM_OVERLAP_NONE 0xFFFF
#define PARAM_OVERLAP_DEFAULT 0
CMachineParameter const paramOverlap =
{
	pt_word,                    // Parameter data type
	"Overlap",                  // Parameter name as its shown in the parameter window
	"FFT windowing overlap",    // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	500,                        // Maximum value
	PARAM_OVERLAP_NONE,         // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_OVERLAP_DEFAULT       // the default slider value
};

#define PARAM_OVERLAPCURVE_IDX 7
#define PARAM_OVERLAPCURVE_NONE 0xFF
#define PARAM_OVERLAPCURVE_DEFAULT 0
CMachineParameter const paramOverlapCurve =
{
	pt_byte,                        // Parameter data type
	"OverlapCurve",                 // Parameter name as its shown in the parameter window
	"FFT windowing overlap curve",  // Parameter description as its shown in the pattern view's statusbar
	0,                              // Minimum value
	OVERLAPCURVE_LINEAR,            // Maximum value
	PARAM_OVERLAPCURVE_NONE,        // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                      // Parameter options, MPF_STATE makes it appears as a slider
	PARAM_OVERLAPCURVE_DEFAULT      // the default slider value
};

#define PARAM_DEBUG_IDX 8
CMachineParameter const paramDebug =
{
	pt_byte,                    // Parameter data type
	"Debug",                    // Parameter name as its shown in the parameter window
	"Debug",                    // Parameter description as its shown in the pattern view's statusbar
	0,                          // Minimum value
	0xFE,                       // Maximum value
	0xFF,                       // Novalue, this value means "nothing happened" in the mi::Tick procedure
	MPF_STATE,                  // Parameter options, MPF_STATE makes it appears as a slider
	0                           // the default slider value
};

CMachineParameter const *pParameters[] = {
	&paramFloor,
	&param2,
	&paramBinSize,
	&paramWindowingFunction,
	&paramAlpha,
	&paramBeta,
	&paramOverlap,
	&paramOverlapCurve,
	&paramDebug
};

CMachineAttribute const *pAttributes[] = { NULL };

#pragma pack(1)                        

class gvals
{
public:
	byte paramFloor;
	word param2;
	byte paramBinSize;
	byte paramWindowingFunction;
	word paramAlpha;
	word paramBeta;
	word paramOverlap;
	byte paramOverlapCurve;
	byte paramDebug;
};

#pragma pack()

CMachineInfo const MacInfo = 
{
	MT_EFFECT,                  // Machine type
	MI_VERSION,                 // Machine interface version
	MIF_DOES_INPUT_MIXING,		// Machine flags
	0,                          // min tracks
	0,                          // max tracks
	PARAM_DEBUG_IDX+1,          // numGlobalParameters
	0,                          // numTrackParameters
	pParameters,				// pointer to parameter stuff
	0,                          // numAttributes
	pAttributes,				// pointer to attribute stuff
	MACHINE_FULL_NAME,          // Full Machine Name
	MACHINE_NAME,               // Short name
	MACHINE_AUTHOR,				// Author name
	"&About..."                 // Right click menu commands
};

class mi;

////////////////////////////////////////////////////////////////////////
// The "miex" class is needed for handling multiple inputs
////////////////////////////////////////////////////////////////////////

class miex : public CMDKMachineInterfaceEx
{
public:
	//Only comment out AddInput, DeleteInput and RenameInput if you also intend to use Input
	//virtual void AddInput(char const *macname, bool stereo) { }
	//virtual void DeleteInput(char const *macename) { }
	//virtual void RenameInput(char const *macoldname, char const *macnewname) { }
	//virtual bool HandleInput(int index, int amp, int pan) { return false; }
	//virtual void SetInputChannels(char const *macname, bool stereo) {}
	//virtual void Input(float *psamples, int numsamples, float amp) {}
	mi *pmi; // pointer to 'parent'
};


////////////////////////////////////////////////////////////////////////
// The main machine interface class
////////////////////////////////////////////////////////////////////////

class mi : public CMDKMachineInterface
{
public:
	mi();
	virtual ~mi();
	virtual void Tick();
	virtual void Stop() {}
	virtual void MDKInit(CMachineDataInput * const pi);
	virtual bool MDKWork(float *psamples, int numsamples, int const mode) { return false; };
	virtual bool MDKWorkStereo(float *psamples, int numsamples, int const mode);
	virtual void MDKSave(CMachineDataOutput * const po) {};
	virtual char const *DescribeValue(int const param, int const value);
	virtual CMDKMachineInterfaceEx *GetEx() { return &ex; }
	virtual void Command(const int i);
	virtual void OutputModeChanged(bool stereo) {}

	void Process(CArray& x);
	void ReallocBuffers();
	void CalculateWindowFunction();
	void CalculateOverlapFunction();
	void ResetOffsets();
	int GetFirstCountDown();
	int GetOverlapSize();
	int GetNextCountDown();
	void CreateOverlapData(FArray& overlapFunction, int length, int overlapType);
public:
	miex ex;
	gvals gval;

	char *debug;
	std::ofstream *debugfile;

	byte pFloorEnabled, p2Enabled;
	float pFloor, p2;
	int pBinSize, pWindowingFunction, pOverlapCurve;
	float pAlpha, pBeta, pOverlap;

	int overlapSize;
	int data1Offset, data2Offset, output1Offset, output2Offset;
	int data1CountDown, data2CountDown;
	CArray data1;
	CArray data2;
	FArray window;
	FArray output1;
	FArray output2;
	FArray overlapCurve;
};

DLL_EXPORTS

mi::mi()
{
	GlobalVals = &gval;
	ex.pmi = this;

	debugfile = new std::ofstream("c:\\temp\\fourier_debug.txt", std::ios::out);

	pFloorEnabled = false;
	p2Enabled = false;
	pBinSize = parameterValueToInt(PARAM_BINSIZE_IDX, PARAM_BINSIZE_DEFAULT);
	pWindowingFunction = parameterValueToInt(PARAM_WINDOWFN_IDX, PARAM_WINDOWFN_DEFAULT);
	pAlpha = parameterValueToInt(PARAM_ALPHA_IDX, PARAM_ALPHA_DEFAULT);
	pBeta = parameterValueToInt(PARAM_BETA_IDX, PARAM_BETA_DEFAULT);
	
	overlapSize = 0;

	debug = new char[16];
	sprintf(debug, "");

	ReallocBuffers(); // also sets outputOffset = -1	
	CalculateWindowFunction(); // also sets data1Offset = 0 

}

mi::~mi()
{
	// Delete allocated memory etc.
	if (debug)
		delete debug;

	if (debugfile)
	{
		debugfile->close();
		delete debugfile;
	}
}

void mi::MDKInit(CMachineDataInput * const pi)
{
	SetOutputMode(true);
	pCB->SetnumOutputChannels(pCB->GetThisMachine(), 2);
}

int _numsamples;

float parameterValueToFloat(int param, word value)
{
	switch (param)
	{
	case PARAM_FLOOR_IDX:
		return (value - 1) / 253.f;
	case PARAM2_IDX:
		return (value - 1);
	case PARAM_ALPHA_IDX:
		return value / 1000.f;
	case PARAM_BETA_IDX:
		return value / 1000.f;
	case PARAM_OVERLAP_IDX:
		return value / 10.f;
	default:
		return (float)value;
	}
}
int parameterValueToInt(int param, word value)
{
	switch (param)
	{
	case PARAM_BINSIZE_IDX:
		return 1 << value;
	default:
		return value;
	}
}

void mi::Tick()
{
	if (gval.paramFloor != PARAM_FLOOR_NONE)
	{
		pFloorEnabled = gval.paramFloor != 0;
		pFloor = parameterValueToFloat(PARAM_FLOOR_IDX, gval.paramFloor);
	}
	if (gval.param2 != PARAM2_NONE)
	{
		p2Enabled = gval.param2 != 0;
		p2 = parameterValueToFloat(PARAM2_IDX, gval.param2);
	}
	if (gval.paramBinSize != PARAM_BINSIZE_NONE)
	{
		pBinSize = parameterValueToInt(PARAM_BINSIZE_IDX, gval.paramBinSize);
	}
	if (gval.paramWindowingFunction != PARAM_WINDOWFN_NONE)
	{
		pWindowingFunction = parameterValueToInt(PARAM_WINDOWFN_IDX, gval.paramWindowingFunction);
	}
	if (gval.paramAlpha != PARAM_ALPHA_NONE)
	{
		pAlpha = parameterValueToFloat(PARAM_ALPHA_IDX, gval.paramAlpha);
	}
	if (gval.paramBeta != PARAM_BETA_NONE)
	{
		pBeta = parameterValueToFloat(PARAM_BETA_IDX, gval.paramBeta);
	}
	if (gval.paramOverlap != PARAM_OVERLAP_NONE)
	{
		pOverlap = parameterValueToFloat(PARAM_OVERLAP_IDX, gval.paramOverlap);
	}
	if (gval.paramOverlapCurve != PARAM_OVERLAPCURVE_NONE)
	{
		pOverlapCurve = parameterValueToInt(PARAM_OVERLAPCURVE_IDX, gval.paramOverlapCurve);
	}

	// Finally
	if (gval.paramOverlap != PARAM_OVERLAP_NONE || gval.paramOverlapCurve != PARAM_OVERLAPCURVE_NONE || gval.paramBinSize != PARAM_BINSIZE_NONE)
	{
		overlapSize = GetOverlapSize();
		if (overlapSize > 0)
			CalculateOverlapFunction();
	}

	if (gval.paramBinSize != PARAM_BINSIZE_NONE)
	{
		ReallocBuffers();
	}

	if (gval.paramBinSize != PARAM_BINSIZE_NONE || gval.paramWindowingFunction != PARAM_WINDOWFN_NONE || gval.paramAlpha != PARAM_ALPHA_NONE || gval.paramBeta != PARAM_BETA_NONE)
	{
		CalculateWindowFunction();
	}
}

int debugT = 0;

bool mi::MDKWorkStereo(float *psamples, int numsamples, const int mode)
{
	// WM_READWRITE when there is input
	// WM_WRITE when there is no input
	
	// If input stops, reset fft buffer
	if (mode == WM_WRITE)
	{
		// Maybe offsets should not be reset just because we get a mode != WM_READWRITE
		ResetOffsets();
	}

	float w, s;
	char txt[64];

	//if (debugfile && mode == 3)
	//{
	//	sprintf(txt, "W %d %.1f (D1 %d %d) (D2 %d %d)\n",  numsamples, pOverlap, data1CountDown, data1Offset, data2CountDown, data2Offset);
	//	*debugfile << txt;
	//}
	if (mode == WM_READWRITE || (mode == WM_WRITE && (output1Offset >= 0 || output2Offset >= 0)))
	{
		for (int i = 0; i < numsamples; ++i)
		{
			if (mode == WM_READWRITE)
			{
				debugT++;
				if (data1CountDown == 0)
				{
					w = window[data1Offset];
					s = 0.5 * (psamples[i * 2] + psamples[(i * 2) + 1]);
					data1[data1Offset] = w * s / 32768.f;
					data1Offset++;
					if (data1Offset == pBinSize)
					{
						//if (debugfile) 
						//{
						//	sprintf(txt, "D1 %d %d\n", debugT - pBinSize, pBinSize);
						//	*debugfile << txt;
						//}

						data1Offset = 0;
						if (overlapSize > 0)
						{
							data1CountDown = GetNextCountDown();
						}

						// forward fft
						fft_fast(data1);

						// process data
						Process(data1);

						// inverse fft
						ifft_fast(data1);

						output1Offset = 0;

						for (int j = 0; j < pBinSize; ++j)
						{
							w = 1.f - window[data1Offset];
							output1[j] = data1[j].real() *32768.f;
						}
					}
				}
				else
				{
					data1CountDown--;
				}

				if (overlapSize > 0)
				{
					if (data2CountDown == 0)
					{
						w = window[data2Offset];
						s = 0.5 * (psamples[i * 2] + psamples[(i * 2) + 1]);
						data2[data2Offset] = w * s / 32768.f;
						data2Offset++;
						if (data2Offset == pBinSize)
						{
							//if (debugfile)
							//{
							//	sprintf(txt, "D2 %d %d\n", debugT - pBinSize, pBinSize);
							//	*debugfile << txt;
							//}

							data2Offset = 0;
							data2CountDown = GetNextCountDown();

							// forward fft
							fft_fast(data2);

							// process data
							Process(data2);

							// inverse fft
							ifft_fast(data2);

							output2Offset = 0;

							for (int j = 0; j < pBinSize; ++j)
							{
								output2[j] = data2[j].real() *32768.f;
							}
						}
					}
					else
					{
						data2CountDown--;
					}
				}
			}

			float k = 0.f;
			bool merge = output1Offset >= 0 && output2Offset >= 0 && overlapSize > 0;
			float kd = merge ? output2Offset < output1Offset ? (1.f - overlapCurve[output2Offset]) : overlapCurve[output1Offset] : 0.f;

			//sprintf(txt, "M %d %d %d %d\n",  overlapSize, output1Offset, output2Offset, merge);
			//*debugfile << txt;

			if (output1Offset >= 0)
			{
				k += merge ? output1[output1Offset] * kd : output1[output1Offset];
				output1Offset++;
				if (output1Offset == pBinSize) {
					output1Offset = -1;
				}
			}
			if (output2Offset >= 0 && overlapSize > 0)
			{
				k += merge ? output2[output2Offset] * (1.f - kd) : output2[output2Offset];
				output2Offset++;
				if (output2Offset == pBinSize) {
					output2Offset = -1;
				}
			}

			psamples[i * 2] = k;
			psamples[(i * 2) + 1] = k;
		}
		return true;
	}

	return false;
}

void mi::Command(int const i)
{
	char message[512];
	char title[32];
	switch (i)
	{
	case 0:
		sprintf(message, "%s v%s\n(c) %s 2019\n\nBlah blah blah", MACHINE_NAME, MACHINE_VERSION, MACHINE_AUTHOR);
		sprintf(title, "About %s", MACHINE_NAME);
		MessageBox(NULL, message, title, MB_OK|MB_SYSTEMMODAL);
		break;
	default:
		break;
	}
}

char const *mi::DescribeValue(const int param, const int value)
{
	static char txt[16];
	switch (param)
	{
	case PARAM_FLOOR_IDX:
		sprintf(txt, value == 0 ? "off" : "%.1f %%", parameterValueToFloat(PARAM_FLOOR_IDX, value)*100);
		return txt;
	case PARAM2_IDX:
		sprintf(txt, value == 0 ? "off" : "%.0f", parameterValueToFloat(PARAM2_IDX, value));
		return txt;
	case PARAM_BINSIZE_IDX:
		sprintf(txt, "%d", parameterValueToInt(PARAM_BINSIZE_IDX, value));
		return txt;
	case PARAM_WINDOWFN_IDX:
		switch (value)
		{
		case WINDOWFN_SQUARE:
			sprintf(txt, "Square"); break;
		case WINDOWFN_HAMMING:
			sprintf(txt, "Hamming"); break;
		case WINDOWFN_HANN:
			sprintf(txt, "Hann"); break;
		case WINDOWFN_GAUSS:
			sprintf(txt, "Gauss"); break;
		case WINDOWFN_FLAT_TOP:
			sprintf(txt, "FlatTop"); break;
		case WINDOWFN_BLACKMAN:
			sprintf(txt, "Blackman"); break;
		case WINDOWFN_BLACKMAN_HARRIS:
			sprintf(txt, "BlackHarris"); break;
		case WINDOWFN_BLACKMAN_NUTTALL:
			sprintf(txt, "BlackNuttall"); break;
		case WINDOWFN_NUTTALL:
			sprintf(txt, "Nuttall"); break;
		case WINDOWFN_KAISER:
			sprintf(txt, "Kaiser"); break;
		case WINDOWFN_KAISER_BESSEL:
			sprintf(txt, "KaiserBessel"); break;
		case WINDOWFN_TRAPEZOID:
			sprintf(txt, "Trapezoid"); break;
		case WINDOWFN_SINC:
			sprintf(txt, "Sinc"); break;
		case WINDOWFN_SINE:
			sprintf(txt, "Sine"); break;
		default:
			return NULL;
		}
		return txt;
	case PARAM_ALPHA_IDX:
		sprintf(txt, "%.3f", parameterValueToFloat(PARAM_ALPHA_IDX, value));
		return txt;
	case PARAM_BETA_IDX:
		sprintf(txt, "%.3f", parameterValueToFloat(PARAM_BETA_IDX, value));
		return txt;
	case PARAM_OVERLAP_IDX:
		sprintf(txt, value == 0 ? "none" : "%.1f %%", parameterValueToFloat(PARAM_OVERLAP_IDX, value));
		return txt;
	case PARAM_OVERLAPCURVE_IDX:
		switch (value)
		{
		case OVERLAPCURVE_EXPONENTIAL:
			sprintf(txt, "Exp"); break;
		case OVERLAPCURVE_LINEAR:
			sprintf(txt, "Linear"); break;
		default:
			return NULL;
		}
		return txt;
	case PARAM_DEBUG_IDX:
		return debug;
	default:
		return NULL;
	}
}

void mi::ReallocBuffers()
{
	data1 = CArray(pBinSize);
	data2 = CArray(pBinSize);
	
	ResetOffsets();

	output1 = FArray(pBinSize);
	output2 = FArray(pBinSize);
	output1Offset = -1;
	output2Offset = -1;
}

void mi::CalculateWindowFunction()
{
	create_window_data(window, pBinSize, pWindowingFunction, pAlpha, pBeta, false);

	ResetOffsets();
}

void mi::CalculateOverlapFunction()
{
	CreateOverlapData(overlapCurve, overlapSize, pOverlapCurve);

	ResetOffsets();
	output1Offset = -1;
	output2Offset = -1;
}

void mi::ResetOffsets()
{
	data1Offset = 0;
	data1CountDown = 0;
	if (overlapSize > 0)
	{
		data2CountDown = GetFirstCountDown();
		data2Offset = 0;
	}
}

int mi::GetFirstCountDown()
{
	return pBinSize - (int)((pBinSize * pOverlap) / 100.f);
}

int mi::GetOverlapSize()
{
	return (int)((pBinSize * pOverlap) / 100.f);
}

int mi::GetNextCountDown()
{
	return pBinSize - 2*GetOverlapSize();
}

void mi::Process(CArray& x)
{
	unsigned int length = x.size();

	float mx = -1.;
	for (unsigned int a = 0; a < length; ++a)
	{
		float fb = fabs(x[a].real());
		if (fb > mx) mx = fb;
		// Lowpass
		//if (a > 20)
		//{
		//	x[a].real(0.f);
		//	x[a].imag(0.f);
		//}
		
		// Floor
		if (pFloorEnabled)
		{
			if (fabs(x[a].real()) < pFloor)
			{
				x[a].real(0.f);
				x[a].imag(0.f);
			}
		}
	}
	sprintf(debug, "%.3f", mx);
}

////////////////////////////////////////////////////////////////////////
// FFT
////////////////////////////////////////////////////////////////////////


// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

// Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
// Better optimized but less intuitive
void fft_fast(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = M_PI / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}

// inverse fft (in-place)
void ifft(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}


void ifft_fast(CArray& x)
{
	// conjugate the complex numbers
	x = x.apply(std::conj);

	// forward fft
	fft_fast(x);

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}

double bessel(double x)
{
	double Sum = 0.0, XtoIpower;
	int i, j, Factorial;
	for (i = 1; i<10; i++)
	{
		XtoIpower = pow(x / 2.0, (double)i);
		Factorial = 1;
		for (j = 1; j <= i; j++)Factorial *= j;
		Sum += pow(XtoIpower / (double)Factorial, 2.0);
	}
	return(1.0 + Sum);
}

double sinc(double x)
{
	if (x > -1.0E-5 && x < 1.0E-5)return(1.0);
	return(sin(x) / x);
}

// See: http://www.iowahills.com/Example%20Code/WindowedFIRFilterWebCode.txt
// These are the various windows definitions. These windows can be used for either
// FIR filter design or with an FFT for spectral analysis.
// Sourced verbatim from: ~MyDocs\Code\Common\FFTFunctions.cpp
// For definitions, see this article:  http://en.wikipedia.org/wiki/Window_function

// This function has 6 inputs
// Data is the array, of length N, containing the data to to be windowed. 
// This data is either a FIR filter sinc pulse, or the data to be analyzed by an fft.

// WindowType is an enum defined in the header file, which is at the bottom of this file.
// e.g. wtKAISER, wtSINC, wtHANNING, wtHAMMING, wtBLACKMAN, ...

// Alpha sets the width of the flat top.
// Windows such as the Tukey and Trapezoid are defined to have a variably wide flat top.
// As can be seen by its definition, the Tukey is just a Hanning window with a flat top.
// Alpha can be used to give any of these windows a partial flat top, except the Flattop and Kaiser.
// Alpha = 0 gives the original window. (i.e. no flat top)
// To generate a Tukey window, use a Hanning with 0 < Alpha < 1
// To generate a Bartlett window (triangular), use a Trapezoid window with Alpha = 0.
// Alpha = 1 generates a rectangular window in all cases. (except the Flattop and Kaiser)


// Beta is used with the Kaiser, Sinc, and Sine windows only.
// These three windows are primarily used for FIR filter design, not spectral analysis.
// In FIR filter design, Beta controls the filter's transition bandwidth and the sidelobe levels.
// The code ignores Beta except in the Kaiser, Sinc, and Sine window cases.

// UnityGain controls whether the gain of these windows is set to unity.
// Only the Flattop window has unity gain by design. The Hanning window, for example, has a gain
// of 1/2.  UnityGain = true will set the gain of all these windows to 1.
// Then, when the window is applied to a signal, the signal's energy content is preserved.
// Don't use this with FIR filter design however. Since most of the enegy in an FIR sinc pulse
// is in the middle of the window, the window needs a peak amplitude of one, not unity gain.
// Setting UnityGain = true will simply cause the resulting FIR filter to have excess gain.

// If using these windows for FIR filters, start with the Kaiser, Sinc, or Sine windows and
// adjust Beta for the desired transition BW and sidelobe levels (set Alpha = 0).
// While the FlatTop is an excellent window for spectral analysis, don't use it for FIR filter design.
// It has a peak amplitude of ~ 4.7 which causes the resulting FIR filter to have about this much gain.
// It works poorly for FIR filters even if you adjust its peak amplitude.
// The Trapezoid also works poorly for FIR filter design.

// If using these windows with an fft for spectral analysis, start with the Hanning, Gauss, or Flattop.
// When choosing a window for spectral analysis, you must trade off between resolution and amplitude accuracy.
// The Hanning has the best resolution while the Flatop has the best amplitude accuracy.
// The Gauss is midway between these two for both accuracy and resolution.
// These three were the only windows available in the HP 89410A Vector Signal Analyzer. Which is to say,
// unless you have specific windowing requirements, use one of these 3 for general purpose signal analysis.
// Set UnityGain = true when using any of these windows for spectral analysis to preserve the signal's enegy level.
void create_window_data(FArray& window, int length, int windowType, double alpha, double beta, bool unityGain)
{
	int j, M, TopWidth;
	double dM;

	if (windowType != WINDOWFN_SQUARE)
	{
		if (windowType == WINDOWFN_KAISER_BESSEL || windowType == WINDOWFN_FLAT_TOP) alpha = 0.0;

		if (alpha < 0.0)alpha = 0.0;
		if (alpha > 1.0)alpha = 1.0;

		if (beta < 0.0)beta = 0.0;
		if (beta > 10.0)beta = 10.0;

		window = FArray(length + 2);

		TopWidth = (int)(alpha * (double)length);
		if (TopWidth % 2 != 0)TopWidth++;
		if (TopWidth > length)TopWidth = length;
		M = length - TopWidth;
		dM = M + 1;


		// Calculate the window for length/2 points, then fold the window over (at the bottom).
		// TopWidth points will be set to 1.
		if (windowType == WINDOWFN_KAISER)
		{
			double Arg;
			for (j = 0; j < M; j++)
			{
				Arg = beta * sqrt(1.0 - pow(((double)(2 * j + 2) - dM) / dM, 2.0));
				window[j] = bessel(Arg) / bessel(beta);
			}
		}

		else if (windowType == WINDOWFN_SINC)  // Lanczos
		{
			for (j = 0; j < M; j++)window[j] = sinc((double)(2 * j + 1 - M) / dM * M_PI);
			for (j = 0; j < M; j++)window[j] = pow(window[j], beta);
		}

		else if (windowType == WINDOWFN_SINE)  // Hanning if beta = 2
		{
			for (j = 0; j < M / 2; j++)window[j] = sin((double)(j + 1) * M_PI / dM);
			for (j = 0; j < M / 2; j++)window[j] = pow(window[j], beta);
		}

		else if (windowType == WINDOWFN_HANN)
		{
			for (j = 0; j < M / 2; j++)window[j] = 0.5 - 0.5 * cos((double)(j + 1) * M_2PI / dM);
		}

		else if (windowType == WINDOWFN_HAMMING)
		{
			for (j = 0; j < M / 2; j++)
				window[j] = 0.54 - 0.46 * cos((double)(j + 1) * M_2PI / dM);
		}

		else if (windowType == WINDOWFN_BLACKMAN)
		{
			for (j = 0; j < M / 2; j++)
			{
				window[j] = 0.42
					- 0.50 * cos((double)(j + 1) * M_2PI / dM)
					+ 0.08 * cos((double)(j + 1) * M_2PI * 2.0 / dM);
			}
		}


		// See: http://www.bth.se/fou/forskinfo.nsf/0/130c0940c5e7ffcdc1256f7f0065ac60/$file/ICOTA_2004_ttr_icl_mdh.pdf
		else if (windowType == WINDOWFN_FLAT_TOP)
		{
			for (j = 0; j <= M / 2; j++)
			{
				window[j] = 1.0
					- 1.93293488969227 * cos((double)(j + 1) * M_2PI / dM)
					+ 1.28349769674027 * cos((double)(j + 1) * M_2PI * 2.0 / dM)
					- 0.38130801681619 * cos((double)(j + 1) * M_2PI * 3.0 / dM)
					+ 0.02929730258511 * cos((double)(j + 1) * M_2PI * 4.0 / dM);
			}
		}


		else if (windowType == WINDOWFN_BLACKMAN_HARRIS)
		{
			for (j = 0; j < M / 2; j++)
			{
				window[j] = 0.35875
					- 0.48829 * cos((double)(j + 1) * M_2PI / dM)
					+ 0.14128 * cos((double)(j + 1) * M_2PI * 2.0 / dM)
					- 0.01168 * cos((double)(j + 1) * M_2PI * 3.0 / dM);
			}
		}

		else if (windowType == WINDOWFN_BLACKMAN_NUTTALL)
		{
			for (j = 0; j < M / 2; j++)
			{
				window[j] = 0.3535819
					- 0.4891775 * cos((double)(j + 1) * M_2PI / dM)
					+ 0.1365995 * cos((double)(j + 1) * M_2PI * 2.0 / dM)
					- 0.0106411 * cos((double)(j + 1) * M_2PI * 3.0 / dM);
			}
		}

		else if (windowType == WINDOWFN_NUTTALL)
		{
			for (j = 0; j < M / 2; j++)
			{
				window[j] = 0.355768
					- 0.487396 * cos((double)(j + 1) * M_2PI / dM)
					+ 0.144232 * cos((double)(j + 1) * M_2PI * 2.0 / dM)
					- 0.012604 * cos((double)(j + 1) * M_2PI * 3.0 / dM);
			}
		}

		else if (windowType == WINDOWFN_KAISER_BESSEL)
		{
			for (j = 0; j <= M / 2; j++)
			{
				window[j] = 0.402
					- 0.498 * cos(M_2PI * (double)(j + 1) / dM)
					+ 0.098 * cos(2.0 * M_2PI * (double)(j + 1) / dM)
					+ 0.001 * cos(3.0 * M_2PI * (double)(j + 1) / dM);
			}
		}

		else if (windowType == WINDOWFN_TRAPEZOID) // Rectangle for alpha = 1  Triangle for alpha = 0
		{
			int K = M / 2;
			if (M % 2)K++;
			for (j = 0; j < K; j++)window[j] = (double)(j + 1) / (double)K;
		}


		// This definition is from http://en.wikipedia.org/wiki/Window_function (Gauss Generalized normal window)
		// We set their p = 2, and use alpha in the numerator, instead of Sigma in the denominator, as most others do.
		// alpha = 2.718 puts the Gauss window response midway between the Hanning and the Flattop (basically what we want).
		// It also gives the same BW as the Gauss window used in the HP 89410A Vector Signal Analyzer.
		// alpha = 1.8 puts it quite close to the Hanning.
		else if (windowType == WINDOWFN_GAUSS)
		{
			for (j = 0; j < M / 2; j++)
			{
				window[j] = ((double)(j + 1) - dM / 2.0) / (dM / 2.0) * 2.7183;
				window[j] *= window[j];
				window[j] = exp(-window[j]);
			}
		}

		else // Error.
		{
			windowType == WINDOWFN_SQUARE;
		}
	}

	if (windowType != WINDOWFN_SQUARE)
	{
		// Fold the coefficients over.
		for (j = 0; j<M / 2; j++)
			window[length - j - 1] = window[j];

		// This is the flat top if alpha > 0. Cannot be applied to a Kaiser or Flat Top.
		if (windowType != WINDOWFN_KAISER &&  windowType != WINDOWFN_FLAT_TOP)
		{
			for (j = M / 2; j<length - M / 2; j++)window[j] = 1.0;
		}

		// This will set the gain of the window to 1. Only the Flattop window has unity gain by design. 
		if (unityGain)
		{
			double Sum = 0.0;
			for (j = 0; j<length; j++)Sum += window[j];
			Sum /= (double)length;
			if (Sum != 0.0)for (j = 0; j<length; j++)window[j] /= Sum;
		}
	}
	else
	{
		window = FArray(length);
		for (j = 0; j < length; j++)
		{
			window[j] = 1.f;
		}
	}
}


//Start with taking the iFFT of S(1) and S(2).Let's call the results m(1) and m(2), since they are modified. They are back in the time domain and you wish to recombine them into one signal again. Ideally, they would line up perfectly and you could just truncate one and keep the other. This won't happen, so the trick is to fade the latter into the former over the overlap region.At the beginning of the region the weighting of m(1) should be one and the weighting of m(2) should be zero.At the end of the region, the reverse should be true.The simplest solution is a linear transition, but that may not be the best.
//
//For a 50 % overlap, there is a nice property of a VonHann window in that the appropriate weighting is implicitely correct in the transition region and you can simply sum the two modified signals over the overlap region.For other overlap sizes the same function can still be used, but it needs to be resized appropriately.The window function can be applied in the time domain before you take the initial FFT, or after the iFFT is done.Of course, which you do greatly affects the bin values, so it depends on the nature of your "touching" on which way is better.The equivalent of a VonHann window can also be done in the frequency domain by calculating a new value for each bin by subtracting the average value of its neighbors and taking half that value.The equation would be Y[k] = ?.25?X[k?1] + .5?X[k]?.25?X[k + 1] where the Xs are the original bin values.Doing this before or after your "touching" is the same as if you were windowing the signal before going into the frequency domain or after coming out.
//
//The less touching you do, the more similar the two modified signals are going to be so the smoother the fading transition.For too much touching, the mismatch between the two modified signals may be too great and the results won't be so good. However, with a fading over the overlap region, the results are always going to be a smooth signal.
//
//The VonHann window is good for spectrogram displays because it spreads the value of a sharp peak(near integer frequencies) to the neighboring bins and reduces the leakage of in between frequencies, so when you do a tone sweep through the frequency range, the display remains somewhat consistent and less frame size dependent.Chances are that this same property will make it desirable for you to do the window function before you do your touching.
//
//Hope this helps,
void mi::CreateOverlapData(FArray& overlapCurve, int length, int overlapType)
{
	overlapCurve = FArray(length);
	if (overlapType == 1)
	{
		// TBD
		for (int j = 0; j < length; ++j)
		{
			overlapCurve[j] = ((float)j) / length;
		}
	}
	else //if (overlapType == 0)
	{
		for (int j = 0; j < length; ++j)
		{
			overlapCurve[j] = ((float)j) / length;
		}
	}
}




//int test_fft()
//{
//	const Complex test[] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
//	CArray data(test, 8);
//
//	// forward fft
//	fft(data);
//
//	std::cout << "fft" << std::endl;
//	for (int i = 0; i < 8; ++i)
//	{
//		std::cout << data[i] << std::endl;
//	}
//
//	// inverse fft
//	ifft(data);
//
//	std::cout << std::endl << "ifft" << std::endl;
//	for (int i = 0; i < 8; ++i)
//	{
//		std::cout << data[i] << std::endl;
//	}
//	return 0;
//}

#pragma optimize ("", on) 
