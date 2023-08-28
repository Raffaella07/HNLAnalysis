#include "../macros/setStyle.C"


TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold
float lumiTextFont   = 42;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "Preliminary";
TString MCText   = "Simulation";
float extraTextFont = 52;  // default is helvetica-italics

// temll_ETt sizes and text offsets with respect to the top frame
// // in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and emll_ETtra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_sqrtS = "";

bool drawLogo      = false;

double SignalWeight(){
	double lumi,n_events;
	lumi = 5.187706291*1e15;
	n_events= 45.16919828*1e15;
	std::string MCpath;

	double crossec_mc;
	crossec_mc = 543100000*1e-12; 

	double weight = lumi/n_events;

	return weight;


}



	void 
CMS_lumi( TPad* pad, int iPeriod, int iPosX,double lumi )
{            
	bool outOfFrame    = false;
	if( iPosX/10==0 ) 
	{
		outOfFrame = true;
	}
	int alignY_=3;
	int alignX_=2;
	if( iPosX/10==0 ) alignX_=1;
	if( iPosX==0    ) alignX_=1;
	if( iPosX==0    ) alignY_=1;
	if( iPosX/10==1 ) alignX_=1;
	if( iPosX/10==2 ) alignX_=2;
	if( iPosX/10==3 ) alignX_=3;
	//if( iPosX == 0  ) relPosX = 0.12;
	int align_ = 10*alignX_ + alignY_;

	float H = pad->GetWh();
	float W = pad->GetWw();
	float l = pad->GetLeftMargin();
	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	float b = pad->GetBottomMargin();
	//  float e = 0.025;

	pad->cd();

	TString lumiText;
	if( iPeriod==1 )
	{
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
	}
	else if ( iPeriod==2 )
	{
		lumiText += lumi_8TeV;
		lumiText += " (8 TeV)";
	}
	else if( iPeriod==3 ) 
	{
		lumiText = lumi_8TeV; 
		lumiText += " (8 TeV)";
		lumiText += " + ";
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
	}
	else if ( iPeriod==4 )
	{
		lumiText += lumi_13TeV;
		lumiText += " (13 TeV)";
	}
	else if ( iPeriod==7 )
	{ 
		if( outOfFrame ) lumiText += "#scale[0.85]{";
		lumiText += lumi_13TeV; 
		lumiText += " (13 TeV)";
		lumiText += " + ";
		lumiText += lumi_8TeV; 
		lumiText += " (8 TeV)";
		lumiText += " + ";
		lumiText += lumi_7TeV;
		lumiText += " (7 TeV)";
		if( outOfFrame) lumiText += "}";
	}
	else if ( iPeriod==12 )
	{
		lumiText += "8 TeV";
	}
	else if ( iPeriod==5 )
	{
		char l[10];
		sprintf(l, "%0.2f", lumi);
		lumiText += std::string(l)+" fb-1 B-parking";//projected to 41.6 fb-1 B-parking" ;
	}

	std::cout << lumiText << endl;

	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);    

	float extraTextSize = extraOverCmsTextSize*cmsTextSize*0.6;

	latex.SetTextFont(42);
	latex.SetTextAlign(31); 
	latex.SetTextSize(lumiTextSize*t*0.6);    
	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

	if( outOfFrame )
	{
		latex.SetTextFont(cmsTextFont);
		latex.SetTextAlign(11); 
		latex.SetTextSize(cmsTextSize*t*0.3);    
		latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
	}

	pad->cd();

	float posX_=0;
	if( iPosX%10<=1 )
	{
		posX_ =   l + relPosX*(1-l-r);
	}
	else if( iPosX%10==2 )
	{
		posX_ =  l + 0.5*(1-l-r);
	}
	else if( iPosX%10==3 )
	{
		posX_ =  1-r - relPosX*(1-l-r);
	}
	float posY_ = 1-t - relPosY*(1-t-b);
	if( !outOfFrame )
	{
		if( drawLogo )
		{
			posX_ =   l + 0.045*(1-l-r)*W/H;
			posY_ = 1-t - 0.045*(1-t-b);
			float xl_0 = posX_;
			float yl_0 = posY_ - 0.15;
			float xl_1 = posX_ + 0.15*H/W;
			float yl_1 = posY_;
			TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
			TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
			pad_logo->Draw();
			pad_logo->cd();
			CMS_logo->Draw("X");
			pad_logo->Modified();
			pad->cd();
		}
		else
		{
			latex.SetTextFont(cmsTextFont);
			latex.SetTextSize(cmsTextSize*t);
			latex.SetTextAlign(align_);
			latex.DrawLatex(posX_, posY_, cmsText);
			if( writeExtraText ) 
			{
				latex.SetTextFont(extraTextFont);
				latex.SetTextAlign(align_);
				latex.SetTextSize(extraTextSize*t);
				latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
			}
		}
	}
	else if( writeExtraText )
	{
		if( iPosX==0) 
		{
			posX_ =   l +  relPosX*(1-l-r);
			posY_ =   1-t+lumiTextOffset*t;
		}
		latex.SetTextFont(extraTextFont);
		latex.SetTextSize(extraTextSize*t);
		latex.SetTextAlign(align_);
		latex.DrawLatex(posX_, posY_, extraText);      
	}
	return;
}

/*void CMS_lumi( TPad* pad, int iPeriod, int iPosX, double lumi )
  {            
  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
  {
  outOfFrame = true;
  }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
//if( iPosX == 0  ) relPosX = 0.12;
int align_ = 10*alignX_ + alignY_;

float H = pad->GetWh();
float W = pad->GetWw();
float l = pad->GetLeftMargin();
float t = pad->GetTopMargin()*0.5;
float r = pad->GetRightMargin();
float b = pad->GetBottomMargin();
//  float e = 0.025;

pad->cd();

TString lumiText;
if( iPeriod==1 )
{
lumiText += lumi_7TeV;
lumiText += " (7 TeV)";
}
else if ( iPeriod==2 )
{
lumiText += lumi_8TeV;
lumiText += " (8 TeV)";
}
else if( iPeriod==3 ) 
{
lumiText = lumi_8TeV; 
lumiText += " (8 TeV)";
lumiText += " + ";
lumiText += lumi_7TeV;
lumiText += " (7 TeV)";
}
else if ( iPeriod==4 )
{
lumiText += lumi_13TeV;
lumiText += " (13 TeV)";
}
else if ( iPeriod==7 )
{ 
if( outOfFrame ) lumiText += "#scale[0.85]{";
lumiText += lumi_13TeV; 
lumiText += " (13 TeV)";
lumiText += " + ";
lumiText += lumi_8TeV; 
lumiText += " (8 TeV)";
lumiText += " + ";
lumiText += lumi_7TeV;
lumiText += " (7 TeV)";
if( outOfFrame) lumiText += "}";
}
else if ( iPeriod==12 )
{
lumiText += "8 TeV";
}
else if ( iPeriod==0 )
{
lumiText += lumi_sqrtS;
}
else if ( iPeriod==5 )
{
	char l[10];
	sprintf(l, "%0.2f", lumi);
	lumiText += l;
	lumiText +=  " fb^{-1}";
}

std::cout << lumiText << std::endl;

TLatex latex;
latex.SetNDC();
latex.SetTextAngle(0);
latex.SetTextColor(kBlack);    

float extraTextSize = extraOverCmsTextSize*cmsTextSize;

latex.SetTextFont(42);
latex.SetTextAlign(31); 
latex.SetTextSize(lumiTextSize*t);    
//	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

outOfFrame = false;
if( outOfFrame )
{
	latex.SetTextFont(cmsTextFont);
	latex.SetTextAlign(11); 
	latex.SetTextSize(cmsTextSize*t);    
	latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
}

pad->cd();

float posX_=0;
if( iPosX%10<=1 )
{
	posX_ =   l + relPosX*(1-l-r);
}
else if( iPosX%10==2 )
{
	posX_ =  l + 0.5*(1-l-r);
}
else if( iPosX%10==3 )
{
	posX_ =  1-r - relPosX*(1-l-r);
}
float posY_ = 1-t - relPosY*(1-t-b);
if( !outOfFrame )
{
	posX_ = r+6*r ;
	posY_ = 1-t - 0.02*(1-t-b);
	if( drawLogo )
	{
		//	posX_ =   r - 0.1*(1-l-r)*W/H;
		//	posY_ = 1-t - 0.045*(1-t-b);
		float xl_0 = posX_;
		float yl_0 = posY_ - 0.15;
		float xl_1 = posX_ + 0.15*H/W;
		float yl_1 = posY_;
		//	TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
		TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
		pad_logo->Draw();
		pad_logo->cd();
		//		CMS_logo->Draw("X");
		pad_logo->Modified();
		pad->cd();
	}
	else
	{
		latex.SetTextFont(cmsTextFont);
		latex.SetTextSize(cmsTextSize*t);
		latex.SetTextAlign(align_);
		latex.DrawLatex(posX_, posY_, cmsText);
		//		latex.DrawLatex(posX_, posY_-posY_*0.49, cmsText);
		if( writeExtraText ) 
		{
			latex.SetTextFont(extraTextFont);
			latex.SetTextAlign(align_);
			latex.SetTextSize(extraTextSize*t);
			//	latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, cmsText);
			//	latex.DrawLatex(posX_, posY_-0.03*posY_, MCText);
			latex.DrawLatex(posX_, posY_-posY_*0.06, extraText);
			latex.SetTextFont(lumiTextFont);
			latex.DrawLatex(posX_*0.95, posY_-posY_*0.03,lumiText+" Bparking");
			outOfFrame = false;

		}
	}
}
else if( writeExtraText )
{
	if( iPosX==0) 
	{
		posX_ =   l +  relPosX*(1-l-r);
		posY_ =   1-t+lumiTextOffset*t;
	}
	latex.SetTextFont(extraTextFont);
	latex.SetTextSize(extraTextSize*t);
	latex.SetTextAlign(align_);
	latex.DrawLatex(posX_, posY_, extraText);      
}
return;
}

*/


void limPlotter(std::string limFile,double mass ,double xlow, double xup, double ylow, double yup, std::string tag, std::string channel,double ef){

	setStyle();
	TH2D * plotter = new TH2D("plotter","plotter",10,xlow,xup,10,ylow,yup);
	TGraph* line_expected = new TGraph(limFile.c_str(),"%lg %lg lg lg lg lg lg ");
	TGraph* line_observed = new TGraph(limFile.c_str(),"%lg %*lg %*lg %*lg %*lg %*lg %lg ");
	TGraphAsymmErrors* g = new TGraphAsymmErrors(limFile.c_str(),"%lg %lg %lg %lg lg lg lg ");
	TGraphAsymmErrors* o = new TGraphAsymmErrors(limFile.c_str(),"%lg %lg %*lg %*lg %lg %lg lg lg ");

	TCanvas* c = new TCanvas("c","c",800,600);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);
	TF1 *scale = new TF1("scale",("y*"+std::to_string(ef)).c_str());
     //   line->Apply(scale);
      //  g->Apply(scale);
      //  o->Apply(scale);
	TLegend* l = new TLegend();
	TLine* thr = new TLine(xlow,1,xup,1);
	thr->SetLineWidth(2);
	c->SetLogy();
	c->SetLogx();
	c->SetBottomMargin(0.25);
//	plotter->GetXaxis()->SetTitle("|V|^{2}");
//	plotter->GetXaxis()->SetTitle("Mass (GeV)");
	plotter->GetXaxis()->SetTitle("#frac{|V_{#mu} V_{e}|^{2}}{|V_{#mu}|^{2}+|V_{e}|^{2}}");
	plotter->GetXaxis()->SetTitleOffset(1.8);
	plotter->GetYaxis()->SetTitle(" 95% C.L limit on #mu");
	line_expected->SetLineColor(kRed);
//	line_expected->SetLineStyle(9);
	line_observed->SetLineColor(kBlack);
	line_observed->SetLineWidth(2);
	line_expected->SetLineWidth(2);
	g->SetFillColor(kGreen-3);
	o->SetFillColor(kOrange);
	l->AddEntry(line_expected,"expected limit","l");
	l->AddEntry(line_observed,"observed limit","l");
	l->AddEntry(g,"68\% expected","f");
	l->AddEntry(o,"95\% expected","f");
	plotter->Draw();
	o->Draw("same3");
	g->Draw("same3");
	line_expected->Draw("sameL");
	line_observed->Draw("sameL");
	thr->Draw("same");	
	l->Draw("same");
	TLatex latex;
	
	gPad->RedrawAxis();
	CMS_lumi( c, 5, 33, 41.6 );
	c->SaveAs(("CombinedLimits/Mass"+std::to_string(mass)+"_"+tag+"_"+channel+".pdf").c_str());




}





