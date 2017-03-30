#include <iostream>
#include <fstream>
#include <list>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <sstream>
#include <cmath>
#include <cairo/cairo.h>
using std::list;
using  std::vector;
#define coeffiecientGradient 0.5

#define widthOfScreen 400
#define heightOfScreen 400
struct dataInput
{
    double x;
    double y;
};
struct dataOutput
{
    double iResult;
};
struct Node
{
    int m_numOutputs;
    //vector<double> m_vec;
    double m_dValResult;
    double m_dValCombined;
    void initializeWeight(int num);
    vector<double> m_vecWeight;
    double transFunction(double dvalCombined);
    Node(){m_dValCombined=0;}
    double m_partialDifferential;

};
void Node::initializeWeight(int numWeight)
{
    for(int i=0;i<numWeight;i++)
    {
        double dTemp=(double)rand()/(double)RAND_MAX;
        m_vecWeight.push_back(dTemp);
    }

}
class Layer
{
public:
    vector<Node> m_vecNodes;
};
class Net
{
public:
    Net(vector<int>& vec)
    {
        for(int i=0;i<vec.size();i++)
        {
            int iNum=vec[i];
            Layer layer;
            if(i==vec.size()-1)
            {

                for(int j=0; j<iNum; j++) //including the bias in each layer except the output layer
                {
                    Node node;
                    node.m_numOutputs=1;
                    node.initializeWeight(0);
                    layer.m_vecNodes.push_back(node);
                }
            }
            else
            {
                for(int j=0; j<=iNum; j++) //including the bias in each layer except the output layer
                {
                    Node node;
                    node.m_numOutputs=vec[i+1];
                    node.initializeWeight(node.m_numOutputs);
                    layer.m_vecNodes.push_back(node);
                }
                Node& ndBias=layer.m_vecNodes.back();
                ndBias.m_dValResult=1;
            }
            m_vecLayer.push_back(layer);
        }
    };
    vector<Layer> m_vecLayer;
    double feedForward(double dValInput1,double dValInput2);
    double logError(double dNetOut,double dOutPut,list<double>&lstError);
    void updateWeight(double targetOut);
};
double getRandomDouble(double fmin,double fmax);
int main()
{
    std::srand(std::time(0));
    std::ifstream ifile("trainingData.txt");
    std::string strLine;
    list<dataInput >lstDataInput;
    list<dataOutput>lstDataOut;
    //stringstream ss;
    int iLogDebug=0;
    while(!ifile.eof())
    {
        std::getline(ifile,strLine);
        std::stringstream ss(strLine);
        std::string strLable;
        ss>>strLable;
        if(strLable.compare("Input:")==0)
        {
            dataInput dataIn;
        //double dvalY;
        ss>>dataIn.x;
        ss>>dataIn.y;
        lstDataInput.push_back(dataIn);
        iLogDebug++;
        }
        else if(strLable.compare("Output:")==0)
        {
            dataOutput dataOut;
            ss>>dataOut.iResult;
            lstDataOut.push_back(dataOut);
        }



    }

    vector<int> vecLayerNumber{2,4,1};
    Net net(vecLayerNumber);
    list<double> lstError;
    auto itOut=lstDataOut.begin();

    int iCountLogDebug=0;
    for(auto it=lstDataInput.begin();it!=lstDataInput.end();it++,itOut++,iCountLogDebug++)
    {
        //it=lstDataInput.begin();
        dataInput input=*it;
        dataOutput output=*itOut;
        double dNetOut=net.feedForward(input.x,input.y);
        net.logError(dNetOut,output.iResult,lstError);
    }

    std::ofstream ofile("ErrorLog.txt");

    for(auto it=lstError.begin();it!=lstError.end();it++)
    {
        ofile<<*it<<std::endl;

    }
    ofile.close();
    //cout << "Hello world!" << endl;
    /////////////////simulation///////////////
    cairo_t *cr;
    cairo_surface_t *surface;
    surface = cairo_image_surface_create (CAIRO_FORMAT_RGB24,widthOfScreen,heightOfScreen);
    cr = cairo_create (surface);
    cairo_new_path (cr);
    cairo_rectangle(cr, 0, 0, widthOfScreen, heightOfScreen);
    cairo_set_source_rgb (cr, 1, 1, 1);
    //cairo_fill(cr);
    cairo_move_to(cr,0,0);
    cairo_line_to(cr,400,400);
    cairo_stroke(cr);
    for(int i=0;i<3000;i++)
    {
        double dCoorX=getRandomDouble(0,widthOfScreen);
        double dCoorY=getRandomDouble(0,heightOfScreen);
        double dNetOut=net.feedForward(dCoorX,dCoorY);
        if(fabs(dNetOut-1)<1e-10)
        {
            cairo_set_source_rgb(cr,0,255,0);
            cairo_rectangle(cr,dCoorX,dCoorY,2,2);
        }
        else
        {
            cairo_set_source_rgb(cr,255,0,0);

            cairo_arc (cr, dCoorX, dCoorY, 2, 0, 2 * M_PI);
        }
        cairo_stroke(cr);

    }
    cairo_surface_write_to_png(surface, "imageSimu.png");

    return 0;
}
double getRandomDouble(double fmin,double fmax)
{
    int iRand=std::rand();
    double temp=fmin+(double)iRand/(double)RAND_MAX*(fmax-fmin);
    return temp;
}
double  Net::feedForward(double dValInput1,double dValInput2)
{
    vector<Layer>& vecLayer=this->m_vecLayer;
    for(int i=0;i<vecLayer.size();i++)
    {
        Layer& layer=vecLayer[i];
        if(i==0)
        {

                Node& node=layer.m_vecNodes[0];
                //node.m_dValCombined=dValInput1;
                node.m_dValResult=dValInput1/widthOfScreen;//node.transFunction(node.m_dValCombined);
                ///////////////////////////////////////
                Node& node2=layer.m_vecNodes[1];
                //node2.m_dValCombined=dValInput2;
                node2.m_dValResult=dValInput2/heightOfScreen;//node.transFunction(node.m_dValCombined);
                ///////////////////////////////////
                Node& nodeBias=layer.m_vecNodes[2];
                nodeBias.m_dValCombined=1;
                nodeBias.m_dValResult=1;


        }
        else if(i==vecLayer.size()-1)
        {
            Layer& layerPrev=vecLayer[i-1];
            for(int j=0; j<layer.m_vecNodes.size(); j++)
            {
                Node& node=layer.m_vecNodes[j];
                node.m_dValCombined=0;
                for(int index=0;index<layerPrev.m_vecNodes.size();index++)
                {
                    Node& ndPrev=layerPrev.m_vecNodes[index];
                    node.m_dValCombined+=ndPrev.m_dValResult*ndPrev.m_vecWeight[j];

                }
                node.m_dValResult=node.transFunction(node.m_dValCombined)>=0.5?1:0;
            }

        }
        else//hidden layer
        {
            Layer& layerPrev=vecLayer[i-1];
            int j=0;
            for(; j<layer.m_vecNodes.size()-1; j++)
            {
                Node& node=layer.m_vecNodes[j];
                node.m_dValCombined=0;
                for(int index=0;index<layerPrev.m_vecNodes.size();index++)
                {
                    Node& ndPrev=layerPrev.m_vecNodes[index];
                    node.m_dValCombined+=ndPrev.m_dValResult*ndPrev.m_vecWeight[j];
                    int iDebug=0;

                }
                node.m_dValResult=node.transFunction(node.m_dValCombined);
            }
            Node& ndBias=layer.m_vecNodes[j];
            ndBias.m_dValResult=1;
        }

    }
    Layer& layLast=vecLayer.back();
    Node& ndLast=layLast.m_vecNodes[0];
    return ndLast.m_dValResult;


}
double Net::logError(double dNetOut,double dOutPut,list<double>& lstError)
{

    double dValDiff=dNetOut-dOutPut;
    lstError.push_back(dValDiff*dValDiff*1/2);
    ////////////////////////////////////update node weight//////////////
   updateWeight(dOutPut);

}
void Net::updateWeight(double dOutPutTarget)
{
    vector<Layer>& vecLayer=m_vecLayer;
     int iLayerNum=vecLayer.size();
    for(int i=iLayerNum-2;i>=0;i--)
    {
        Layer& layer=vecLayer[i];
        Layer& layerNext=vecLayer[i+1];
        int numNodes=layer.m_vecNodes.size();
        int numNodesNextLayer=layerNext.m_vecNodes.size();
        if(i==iLayerNum-2)
        {
            double deltaWeight=0;
            for(int j=0; j<numNodes; j++)
            {
                Node& nd=layer.m_vecNodes[j];
                nd.m_partialDifferential=0;

                for(int indexWeight=0; indexWeight<nd.m_vecWeight.size(); indexWeight++)
                {
                    //////////////delata Weight//////////////////////////
                    Node& nodeOut=layerNext.m_vecNodes[indexWeight];
                    double fval=nodeOut.transFunction(nodeOut.m_dValCombined);
                    double derivative=(1-fval)*fval;
                    deltaWeight=(nodeOut.m_dValResult-dOutPutTarget)*derivative*nd.m_dValResult;
                    if(fabs(deltaWeight)>1e-10)
                    {
                        int iDebugNow=0;

                    }
                    nd.m_vecWeight[indexWeight]-=coeffiecientGradient*deltaWeight;
                    nd.m_partialDifferential+=deltaWeight;
                }


            }

        }
        else//first layer and hidden layer
        {
            double deltaWeight=0;
            for(int j=0; j<numNodes; j++)
            {
                Node& nd=layer.m_vecNodes[j];
                for(int indexWeight=0; indexWeight<nd.m_vecWeight.size(); indexWeight++)
                {
                    //////////////delata Weight//////////////////////////
                    Node& nodeOut=layerNext.m_vecNodes[indexWeight];
                    double fval=nodeOut.transFunction(nodeOut.m_dValCombined);
                    double derivative=(1-fval)*fval;
                    deltaWeight=(nodeOut.m_partialDifferential)*derivative*nd.m_dValResult;
                    if(fabs(deltaWeight)>1e-10)
                    {
                        int iDebugNow=0;

                    }
                    nd.m_vecWeight[indexWeight]-=coeffiecientGradient*deltaWeight;
                }

            }

        }

    }


}
double Node::transFunction(double dValCombined)
{
    double dTemp=exp(-dValCombined);
    //return tanh(dValCombined);
    return 1/(double)(1+dTemp);
}
