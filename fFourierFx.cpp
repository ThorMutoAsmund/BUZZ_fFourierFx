#include "mdk/mdk.h"
#include "math.h"
#include <windows.h>

#pragma optimize ("awy", on) 

#define PROGRAM_NAME "fFourierFx"
#ifdef _DEBUG
#define PROGRAM_VERSION "0.01.0 beta"
#else
#define PROGRAM_VERSION "1.06.0"
#endif
#define PARAM1_NONE 0xFF
#define PARAM1_DEFAULT 0x00
#define PARAM2_NONE 0xFFFF
#define PARAM2_DEFAULT 0x0000


CMachineParameter const param1 =
{
	pt_byte,                        // Parameter data type
	"P1",                        // Parameter name as its shown in the parameter
                                                // window
	"Parameter 1",            // Parameter description as its shown in
                                                //the pattern view's statusbar
	0,                                    // Minimum value
	0xFE,                        // Maximum value
	PARAM1_NONE,                        // Novalue, this value means "nothing
                                                // happened" in the mi::Tick procedure
	MPF_STATE,                        // Parameter options, MPF_STATE makes it 
                                                // appears as a slider
	PARAM1_DEFAULT                        // the default slider value
};

CMachineParameter const param2 =
{
	pt_word,                        // Parameter data type
	"P2",                        // Parameter name as its shown in the parameter
                                                // window
	"Parameter 2",            // Parameter description as its shown in
                                                //the pattern view's statusbar
	0,                                    // Minimum value
	0x8000,                        // Maximum value
	PARAM2_NONE,                        // Novalue, this value means "nothing
                                                // happened" in the mi::Tick procedure
	MPF_STATE,                        // Parameter options, MPF_STATE makes it 
                                                // appears as a slider
	PARAM2_DEFAULT                        // the default slider value
};



CMachineParameter const *pParameters[] = {
	&param1,
	&param2
};

CMachineAttribute const *pAttributes[] = { NULL };


///////////////////////////////////////////////
// PART TWO
///////////////////////////////////////////////

#pragma pack(1)                        

class gvals
{
public:
	byte param1;
	word param2;
};

#pragma pack()

CMachineInfo const MacInfo = 
{
	MT_EFFECT,                        // Machine type
	MI_VERSION,                        // Machine interface version
	MIF_DOES_INPUT_MIXING,								// Machine flags
	0,                                    // min tracks
	0,                                    // max tracks
	2,                                    // numGlobalParameters
	0,                                    // numTrackParameters
	pParameters,            // pointer to parameter stuff
	0,                                    // numAttributes
	pAttributes,            // pointer to attribute stuff
#ifdef _DEBUG
	"Nafai's fFourierFx (debug)",            // Full Machine Name
#else
	"Nafai's fFourierFx",            // Full Machine Name
#endif
	"fFourierFx",                        // Short name
	"Thor Muto Asmund",            // Author name
	"&About..."                        // Right click menu commands
};


class miex;

class mi;

class miex : public CMDKMachineInterfaceEx
{
public:
	virtual void Input(float *psamples, int numsamples, float amp);
	virtual bool HandleInput(int index, int amp, int pan); //{ return false; }
	virtual void SetInputChannels(char const *macname, bool stereo); //{}
/*
	virtual void AddInput(char const *macname, bool stereo);            // called when input is added to a machine
	virtual void DeleteInput(char const *macename);                                    
	virtual void RenameInput(char const *macoldname, char const *macnewname); 
*/
public:
	mi *pmi;
};


class mi : public CMDKMachineInterface
{
public:
	mi();
	virtual ~mi();
	virtual void Tick();
	virtual void MDKInit(CMachineDataInput * const pi);
	virtual bool MDKWork(float *psamples, int numsamples, int const mode);
	virtual bool MDKWorkStereo(float *psamples, int numsamples, int const mode);
	virtual void Command(int const i);
	virtual void MDKSave(CMachineDataOutput * const po);
	virtual char const *DescribeValue(int const param, int const value);
	virtual CMDKMachineInterfaceEx *GetEx() { return &ex; }
	virtual void OutputModeChanged(bool stereo) {}

public:
	miex ex;

public:
	float p1;
	int arraysize;
	float *thearray;
	//float *fadearray;
	//byte previnputsel,inputsel;
	//float crossfade;
	//int curinput;
	//int fadeval, maxfadeval;

	gvals gval;
};



///////////////////////////////////////////////
// PART THREE
///////////////////////////////////////////////

// Miex
/*
void miex::AddInput(char const *macname, bool stereo) { }
void miex::DeleteInput(char const *macename) { }
void miex::RenameInput(char const *macoldname, char const *macnewname) { }
	*/
void miex::SetInputChannels(char const *macname, bool stereo) { }
bool miex::HandleInput(int index, int amp, int pan) { return false; }
void miex::Input(float *psamples, int numsamples, float amp)
{
	if (numsamples > pmi->arraysize)
	{
		if (pmi->thearray)
			delete pmi->thearray;
		pmi->arraysize = numsamples;
		pmi->thearray = new float[2 * pmi->arraysize];
		for (int i = 0; i < numsamples * 2; i++) {
			pmi->thearray[i] = 0.0f;
		}
	}

	for (int i = 0; i < numsamples * 2; i++)
		pmi->thearray[i] = pmi->thearray[i] + amp*(*psamples++);
}

// Mi
mi::mi() 
{
	GlobalVals = &gval;

	// Init local fields
	thearray = NULL;
	p1 = 0.f;
}

mi::~mi()
{
	// Delete allocated memory in local fields etc
	if (thearray)
		delete thearray;
}

void mi::MDKInit(CMachineDataInput * const pi)
{
	SetOutputMode( true ); // No mono sounds
	ex.pmi = this;
	arraysize = 1024;
	thearray = new float[2 * arraysize];
	/*inputsel = INPUTSEL_DEFAULT;
	previnputsel = INPUTSEL_DEFAULT;
	crossfade = docrossfade(CROSSFADE_DEFAULT);
	arraysize = 1024;
	thearray = new float[2*arraysize];
	fadearray = new float[2*arraysize];
	fadeval = 0;
	maxfadeval = 0;*/
}

void mi::MDKSave(CMachineDataOutput * const po) { }

void mi::Tick() {
	if (gval.param1 != PARAM1_NONE)
	{
		p1 = gval.param1 / 254.f;
		//previnputsel = inputsel;
		//inputsel = gval.inputsel;
		//int storemaxfadeval = maxfadeval;
		//maxfadeval = (int)(crossfade/1000 * pMasterInfo->SamplesPerSec);
		//if (storemaxfadeval == 0)
		//	fadeval = maxfadeval;
		//else
		//	fadeval = (int) (((float)storemaxfadeval-fadeval)/storemaxfadeval*maxfadeval);
	}
	if (gval.param2 != PARAM2_NONE)
	{
		//crossfade = docrossfade(gval.crossfade);
	}
}

bool mi::MDKWork(float *psamples, int numsamples, int const mode)
{
	return false;
}

bool mi::MDKWorkStereo(float *psamples, int numsamples, int const mode)
{
	if (mode==WM_WRITE)
		return false;
	if (mode==WM_NOIO)
		return false;
	if (mode==WM_READ)                        // <thru>
		return false;

	for (int i = 0; i < numsamples; i++)
	{
		*psamples++ = thearray[i * 2] * p1;
		*psamples++ = thearray[i * 2 + 1] * p1;
		thearray[i * 2] = 0.0f;
		thearray[i * 2 + 1] = 0.0f;
	}
	return true;
}

void mi::Command(int const i)
{
	char txt[512];
	switch (i)
	{
	case 0:
		sprintf(txt,"%s v%s\n(c) Thor Muto Asmund 2019\n\nBlah blah blah",PROGRAM_NAME,PROGRAM_VERSION);
		MessageBox(NULL,txt,"About fFourierFx",MB_OK|MB_SYSTEMMODAL);
		break;
	default:
		break;
	}
}
char const *mi::DescribeValue(int const param, int const value)
{
	//static char txt[16];
	switch(param)
	{
	//case 0:
	//	if (value == 0)
	//		strcpy(txt,"None");
	//	else if (value == 17)
	//		strcpy(txt,"All");
	//	else
	//		sprintf(txt,"%d", value );
	//	return txt;
	//	break;
	//case 1:
	//	if (value == 0)
	//		strcpy(txt,"Off");
	//	else
	//		sprintf(txt,"%.1f ms", docrossfade((float)value));
	//	return txt;
	//	break;
	default:
		return NULL;
	}
}

#pragma optimize ("", on) 

DLL_EXPORTS
