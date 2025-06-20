#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/ipv4-global-routing-helper.h"

//
#include "ns3/traffic-control-module.h"
#include "ns3/mobility-module.h"
#include <random> 
#include <unordered_set> 

#include <iomanip>


using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("FatTreeFCT");

enum PROTOCOL {
    TcpCubic, 
    TcpNewReno, 
    TcpBbr, 
    TcpDctcp,
    ProtocolCount
};

std::string protocolName;

// Helper function to install a TCP congestion control algorithm
void SetTcpCongestionControl(int protocol)
{
    switch (protocol)
    {
    case PROTOCOL::TcpCubic:
        protocolName = "TcpCubic";
        break;
    case PROTOCOL::TcpNewReno:
        protocolName = "TcpNewReno";
        break;
    case PROTOCOL::TcpBbr:
        protocolName = "TcpBbr";
        break;
    default:
        protocolName = "TcpDctcp";
        break;
    }

    // Setting the protocol.
    Config::SetDefault("ns3::TcpL4Protocol::SocketType", StringValue("ns3::" + protocolName));
    std::cout << std::endl << "Protocol: " << protocolName << std::endl;
}

using Pair = std::pair<int, int>;

bool traceQueue = false;
int queueThreshold = 5;

//

std::ofstream rttFile;
// rttFile("scratch/rtt.txt", std::ios::out);

void RttTracer(Time oldRtt, Time newRtt)
{
    double now = Simulator::Now().GetSeconds();

    // rttFile << std::fixed << std::setprecision(1) << Simulator::Now().GetSeconds() << "\t" << newRtt.GetMilliSeconds() << std::endl;
    rttFile << now << "\t" << newRtt.GetMilliSeconds() << std::endl;
    
    std::cout << "[RTT TRACE] Time: " << now << "s, RTT Changed: "
              << oldRtt.GetMilliSeconds() << " ms → "
              << newRtt.GetMilliSeconds() << " ms" << std::endl;
}

std::ofstream queueFile;
// queueFile("scratch/queue_size.txt", std::ios::out);
void QueueSizeTracer(uint32_t oldValue, uint32_t newValue)
{
    uint32_t diff = (newValue > oldValue) ? (newValue - oldValue) : (oldValue - newValue);
    if (diff < queueThreshold){
        return;
    }

    double now = Simulator::Now().GetSeconds();
    queueFile << now << "\t" << newValue << std::endl;
    std::cout << "[QUEUE SIZE] Time: " << now 
              << "s, Queue Size Changed: " 
              << oldValue << " → " << newValue << " packets" << std::endl;
}

void PeriodicQueueLogger(Ptr<QueueDisc> queue, int index, int queueNum, double interval) {
    double now = Simulator::Now().GetSeconds();
    uint32_t qSize = queue->GetCurrentSize().GetValue();
    queueFile << index << "\t" << queueNum << "\t" << now << "\t" << qSize << std::endl;
    // Schedule next poll
    Simulator::Schedule(Seconds(interval), &PeriodicQueueLogger, queue, index, queueNum, interval);
}
//

std::ofstream cwndFile;
// cwndFile("scratch/cwnd.txt", std::ios::out);
void CwndTracer(uint32_t oldCwnd, uint32_t newCwnd) 
{
    double time = Simulator::Now().GetSeconds(); // Get simulation time

    // cwndFile << std::fixed << std::setprecision(1) << Simulator::Now().GetSeconds() << "\t" << newCwnd << std::endl; // Write time and cwnd to file
    cwndFile << time << "\t" << newCwnd << std::endl; // Write time and cwnd to file
    std::cout << "[DEBUG] Time: " << time << "s, CWND Changed: " << oldCwnd << " -> " << newCwnd << std::endl;
}

void Track_Flow(Ptr<TcpSocketBase> socket)
{
    socket->TraceConnectWithoutContext("CongestionWindow", MakeCallback(&CwndTracer));

    //
    socket->TraceConnectWithoutContext("RTT", MakeCallback(&RttTracer));
}

int GetRandomNumber(int n) 
{
    std::random_device rd;   
    std::mt19937 gen(rd());     
    std::uniform_int_distribution<> distrib(0, n); // Range [0, n]
    return distrib(gen);       
}

std::vector<int> GetRandomNumbers(int n, int x, int q) 
{
    std::random_device rd;     
    std::mt19937 gen(rd());    
    std::uniform_int_distribution<> distrib(0, x); // Range [0, x]

    std::vector<int> result;
    for (int i = 0; i < n; ++i) 
    {
        int num;
        do {
            num = distrib(gen);
        } while (num == q); // Skip 'q' if you still want to avoid it
        result.push_back(num);
    }
    return result;
}

std::vector<Pair> GetRandomPairs(int n, int x) 
{
    if (n <= 0 || x < 1) { // x must be at least 1 to ensure two distinct numbers
        throw std::invalid_argument("n must be > 0 and x must be >= 1.");
    }
    std::random_device rd;    
    std::mt19937 gen(rd());    
    std::uniform_int_distribution<> dist(0, x); // Range [0, x]
    std::vector<Pair> randomPairs;

    for (int i = 0; i < n; ++i) 
    {
        int first = dist(gen); 
        int second;
        do 
        {
            second = dist(gen);
        } while (second == first);

        randomPairs.emplace_back(first, second);
    }
    return randomPairs;
}

// DCTCP ECN:
static uint32_t g_ecnMarkCount = 0;
std::ofstream ecnFile;

// void EcnMarkTracer(Ptr<const QueueDiscItem> item, const char* reason)
// {
//     g_ecnMarkCount++;
//     // Optionally, you can print/log info about the item or reason
// }

void EcnMarkTracer(Ptr<const QueueDiscItem> item, const char* reason) {
    double now = Simulator::Now().GetSeconds();
    static uint32_t ecnMarkCount = 0; // Or make this global if you want to access elsewhere
    ecnMarkCount++;
    ecnFile << now << "\t" << ecnMarkCount << std::endl;
    std::cout << "[ECN MARK] Time: " << now << "s, ECN Mark Count: " << ecnMarkCount << " Reason: " << reason << std::endl;
}

void PeriodicEcnLogger(double interval) {
    static uint32_t lastMarkCount = 0;
    double now = Simulator::Now().GetSeconds();
    uint32_t diff = g_ecnMarkCount - lastMarkCount;
    ecnFile << now << "\t" << g_ecnMarkCount << "\t" << diff << std::endl;
    lastMarkCount = g_ecnMarkCount;
    Simulator::Schedule(Seconds(interval), &PeriodicEcnLogger, interval);
}



//

/*
k                           Fat-tree switch port count (e.g., 4)
p2p_DataRate:               Link bandwidth (e.g., 10Gbps)
                            the maximum bandwidth of the network link ("pipe")
p2p_Delay:                  Link delay (e.g., 10us)
mode:                       Flow behavior:
                            0 = All to last host
                            1 = All to one random host
                            2 = Fully random
num_of_flows:               Number of concurrent flows (e.g., 100)
flowmonitor_start_time:     Simulation start time (e.g., 0)
simulation_stop_time:       Simulation duration (e.g., 15)
app_data_rate:              Flow-level app traffic rate (e.g., 100Mbps)
                            how fast the application tries to send data ("faucet")
app_packet_size:            Application packet size in bytes (e.g., 1000)

*/

void runGit(int k, std::string p2p_DataRate, std::string p2p_Delay, int mode, uint32_t num_of_flows, 
            double flowmonitor_start_time, double simulation_stop_time,
            std::string app_data_rate, int app_packet_size)
{
    NodeContainer coreSwitches, aggSwitches, edgeSwitches, hosts;

    NS_LOG_INFO("Starting the simulation...");
    
    // Computes the number of core, aggregation, edge switches, and hosts based on the input k.
    // These define the size of the Fat-Tree data center topology.
    uint32_t numPods = k;
    uint32_t numCoreSwitches = (k / 2) * (k / 2);
    uint32_t numAggSwitches = (k * k) / 2;
    uint32_t numEdgeSwitches = numAggSwitches;
    uint32_t numHosts = k * k * k / 4;
    uint32_t threshold_k = 2666;   ////// needs to be checked
    /// uint32_t threshold_k = 250;   ////// needs to be checked

    if (protocolName == "TcpDctcp"){
        /*
        DCTCP is designed to react to early signs of congestion, so we want ECN marking to start when the queue is still shallow. 
        This enables DCTCP to keep queue occupancy low and avoid bufferbloat.
        */
        threshold_k = 20;
    }
    else
    {
        Config::SetDefault("ns3::TcpSocketState::EnablePacing",  BooleanValue(false));
        Config::SetDefault("ns3::TcpSocketState::MaxPacingRate", DataRateValue(DataRate("1Gb/s")));
    }

    //
    bool enableSwitchEcn;
    bool enable_swift;
    //
    
    std::vector<NetDeviceContainer> link_vec;
    
    Config::SetDefault("ns3::TcpL4Protocol::RecoveryType", StringValue("ns3::TcpClassicRecovery"));

    //tcp config 
    // Config::SetDefault("ns3::TcpSocket::SegmentSize", UintegerValue(1500));
    // NEW
    Config::SetDefault("ns3::TcpSocket::SegmentSize", UintegerValue(app_packet_size));

    Config::SetDefault("ns3::TcpSocket::DelAckCount", UintegerValue(1)); // for two packets return one ack
    GlobalValue::Bind("ChecksumEnabled", BooleanValue(false));
    
    //RED config

    // Set default parameters for RED queue disc
    Config::SetDefault("ns3::RedQueueDisc::UseEcn", BooleanValue(true));
    // ARED may be used but the queueing delays will increase; it is disabled
    // here because the SIGCOMM paper did not mention it
    Config::SetDefault ("ns3::RedQueueDisc::Gentle", BooleanValue (false));
    Config::SetDefault("ns3::RedQueueDisc::UseHardDrop", BooleanValue(false));

    // Config::SetDefault("ns3::RedQueueDisc::MeanPktSize", UintegerValue(1500));
    // NEW
    Config::SetDefault("ns3::RedQueueDisc::MeanPktSize", UintegerValue(app_packet_size));

    // Why 2666? The If every packet is 1500 bytes (max) and we want 4MB max, it means up to 2666p.
    
    Config::SetDefault("ns3::RedQueueDisc::MaxSize", QueueSizeValue(QueueSize("2666p")));
    // DCTCP tracks instantaneous queue length only; so set QW = 1
    Config::SetDefault("ns3::RedQueueDisc::QW", DoubleValue(1));
    // in the paper they have said the min max should be equal
    Config::SetDefault("ns3::RedQueueDisc::MinTh", DoubleValue(threshold_k));
    Config::SetDefault("ns3::RedQueueDisc::MaxTh", DoubleValue(threshold_k));


    TrafficControlHelper tchRed1;
    // MinTh = 20, MaxTh = 60 recommended in ACM SIGCOMM 2010 DCTCP Paper
    // This yields a target queue depth of 250us at 1 Gb/s
    tchRed1.SetRootQueueDisc("ns3::RedQueueDisc",
                             "LinkBandwidth",
                             StringValue(p2p_DataRate),
                             "LinkDelay",
                             StringValue(p2p_Delay),
                             "MinTh",
                             DoubleValue(threshold_k),
                             "MaxTh",
                             DoubleValue(threshold_k));

    // Nodes
    coreSwitches.Create(numCoreSwitches);
    aggSwitches.Create(numAggSwitches);
    edgeSwitches.Create(numEdgeSwitches);
    hosts.Create(numHosts);

    PointToPointHelper p2p;
    p2p.SetDeviceAttribute("DataRate", StringValue(p2p_DataRate)); 
    p2p.SetChannelAttribute("Delay", StringValue(p2p_Delay));

    // Internet Stack
    InternetStackHelper internet;
    internet.Install(coreSwitches);
    internet.Install(aggSwitches);
    internet.Install(edgeSwitches);
    internet.Install(hosts);

    int link_index = 0;
    // Core to Aggregation Links
    for (uint32_t pod = 0; pod < numPods; ++pod) {
        for (uint32_t i = 0; i < k / 2; ++i) {
            for (uint32_t j = 0; j < k / 2; ++j) {
                uint32_t coreIndex = i * (k / 2) + j;
                uint32_t aggIndex = pod * (k / 2) + i;
                NetDeviceContainer link = p2p.Install(coreSwitches.Get(coreIndex), aggSwitches.Get(aggIndex));
                link_vec.push_back(link);

                // tchRed1.Install(link.Get(1));
                // tchRed1.Install(link.Get(0));

                // QueueDiscContainer q0 = tchRed1.Install(link.Get(0));
                // QueueDiscContainer q1 = tchRed1.Install(link.Get(1));

                // if (traceQueue)
                // {
                //     Ptr<QueueDisc> disc0 = q0.Get(0);
                //     disc0->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));

                //     Ptr<QueueDisc> disc1 = q1.Get(0);
                //     disc1->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));
                // }


                Ptr<QueueDisc> queueDisc1 = tchRed1.Install(link.Get(1)).Get(0);
                Ptr<QueueDisc> queueDisc0 = tchRed1.Install(link.Get(0)).Get(0);

                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc0, 0, link_index, 0.05); // logs every 0.05s
                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc1, 1, link_index, 0.05); // logs every 0.05s

                if (protocolName == "TcpDctcp") {
                    // queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, queueDisc1, 0.05); // log every 0.05s (change interval as needed)
                    queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, 0.05);
                }

                link_index++;
            }
        }
    }

    link_index = 0;
    // Aggregation to Edge Links
    for (uint32_t pod = 0; pod < numPods; ++pod) {
        for (uint32_t i = 0; i < k / 2; ++i) {
            for (uint32_t j = 0; j < k / 2; ++j) {
                uint32_t aggIndex = pod * (k / 2) + i;
                uint32_t edgeIndex = pod * (k / 2) + j;
                NetDeviceContainer link = p2p.Install(aggSwitches.Get(aggIndex), edgeSwitches.Get(edgeIndex));
                link_vec.push_back(link);
                // tchRed1.Install(link.Get(1));
                // tchRed1.Install(link.Get(0));

                // QueueDiscContainer q0 = tchRed1.Install(link.Get(0));
                // QueueDiscContainer q1 = tchRed1.Install(link.Get(1));

                // if (traceQueue)
                // {
                //     Ptr<QueueDisc> disc0 = q0.Get(0);
                //     disc0->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));

                //     Ptr<QueueDisc> disc1 = q1.Get(0);
                //     disc1->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));
                // }


                Ptr<QueueDisc> queueDisc1 = tchRed1.Install(link.Get(1)).Get(0);
                Ptr<QueueDisc> queueDisc0 = tchRed1.Install(link.Get(0)).Get(0);

                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc0, 0, link_index, 0.05); // logs every 0.05s
                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc1, 1, link_index, 0.05); // logs every 0.05s

                if (protocolName == "TcpDctcp") {
                    // queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, queueDisc1, 0.05); // log every 0.05s (change interval as needed)
                    queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, 0.05);
                }

                link_index++;
            }
        }
    }

    link_index = 0;
    // Edge to Host Links
    for (uint32_t pod = 0; pod < numPods; ++pod) {
        for (uint32_t edgeIndex = pod * (k / 2); edgeIndex < (pod + 1) * (k / 2); ++edgeIndex) {
            for (uint32_t i = 0; i < k / 2; ++i) {
                uint32_t hostIndex = edgeIndex * (k / 2) + i;
                NetDeviceContainer link = p2p.Install(edgeSwitches.Get(edgeIndex), hosts.Get(hostIndex));
                link_vec.push_back(link);
                // tchRed1.Install(link.Get(1));
                // tchRed1.Install(link.Get(0));

                // QueueDiscContainer q0 = tchRed1.Install(link.Get(0));
                // QueueDiscContainer q1 = tchRed1.Install(link.Get(1));

                // if (traceQueue)
                // {
                //     Ptr<QueueDisc> disc0 = q0.Get(0);
                //     disc0->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));

                //     Ptr<QueueDisc> disc1 = q1.Get(0);
                //     disc1->TraceConnectWithoutContext("PacketsInQueue", MakeCallback(&QueueSizeTracer));
                // }


                Ptr<QueueDisc> queueDisc1 = tchRed1.Install(link.Get(1)).Get(0);
                Ptr<QueueDisc> queueDisc0 = tchRed1.Install(link.Get(0)).Get(0);

                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc0, 0, link_index, 0.05); // logs every 0.05s
                Simulator::Schedule(Seconds(2.0), &PeriodicQueueLogger, queueDisc1, 1, link_index, 0.05); // logs every 0.05s

                if (protocolName == "TcpDctcp") {
                    // queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, queueDisc1, 0.05); // log every 0.05s (change interval as needed)
                    queueDisc1->TraceConnectWithoutContext("Mark", MakeCallback(&EcnMarkTracer));
                    // Simulator::Schedule(Seconds(2.0), &PeriodicEcnLogger, 0.05);
                }

                link_index++;
            }
        }
    }

    // Assign IP Addresses
    Ipv4AddressHelper ipv4;
    ipv4.SetBase("10.0.0.0", "255.255.255.0");

    for (uint32_t pod = 0; pod < link_vec.size(); ++pod) {    
        Ipv4InterfaceContainer interface = ipv4.Assign(link_vec[pod]);
        ipv4.NewNetwork();
    }
   
    // Enable Routing
    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    MobilityHelper mobility;
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(coreSwitches);
    mobility.Install(aggSwitches);
    mobility.Install(edgeSwitches);
    mobility.Install(hosts);
    

    //---------------------here starts the implementation of the connections between servers------------------
    uint16_t basePort = 50000; // Base port for applications
    std::vector<ApplicationContainer> sinkApps;
    std::vector<ApplicationContainer> sourceApps;

    // All to last host
    if(mode == 0) 
    {
        // int random_window_to_track  = GetRandomNumber(num_of_flows-1); 
        for (uint32_t i = 0; i < num_of_flows; ++i) 
        {
            uint16_t port = basePort + i;
            if ((i % hosts.GetN()) == ((hosts.GetN() - 1) % hosts.GetN()))
            {
                // Skipping the server.
                continue;
            }
            Ptr<Node> sender = hosts.Get(i % hosts.GetN());
            Ptr<Node> receiver = hosts.Get((hosts.GetN() - 1) % hosts.GetN());


            // Install PacketSink on the receiver
            Address sinkAddress(InetSocketAddress(receiver->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal(), port));
            PacketSinkHelper sinkHelper("ns3::TcpSocketFactory", sinkAddress);
            ApplicationContainer sinkApp = sinkHelper.Install(receiver);
            sinkApp.Start(Seconds(1.0));
            sinkApp.Stop(Seconds(simulation_stop_time));
            sinkApps.push_back(sinkApp);

            //  Install on off on the sender
            OnOffHelper onOffHelper("ns3::TcpSocketFactory", sinkAddress);
            onOffHelper.SetAttribute("DataRate", StringValue(app_data_rate)); // Sending rate
            onOffHelper.SetAttribute("PacketSize", UintegerValue(app_packet_size)); // Packet size (bytes)
            ApplicationContainer sourceApp = onOffHelper.Install(sender);
            sourceApp.Start(Seconds(2.0));
            sourceApp.Stop(Seconds(simulation_stop_time - 1));


            if (i == 0) { // Track the first window
                Ptr<Application> app = sourceApp.Get(0); // Get the OnOffApplication instance
                Ptr<OnOffApplication> onOffApp = DynamicCast<OnOffApplication>(app);
                
                if (onOffApp) {
                    Simulator::Schedule(Seconds(2.01), [onOffApp]() {
                        Ptr<Socket> onOffSocket = onOffApp->GetSocket();
                        Ptr<TcpSocketBase> tcpSocketBase = DynamicCast<TcpSocketBase>(onOffSocket);
                        
                        if (tcpSocketBase) {
                            std::cout << "[DEBUG] Attaching cwnd tracer to connection 0 (Sender Node " << onOffApp->GetNode()->GetId() << ")" << std::endl;
                            Track_Flow(tcpSocketBase);  
                        } else {
                            std::cerr << "[ERROR] Could not cast OnOff socket to TcpSocketBase!" << std::endl;
                        }
                    });
                } else {
                    std::cerr << "[ERROR] Could not retrieve OnOffApplication instance!" << std::endl;
                }
            }
        }
    }
    
    // All to one random host
    else if(mode == 1)
    {  
        // Random a server
        int reciever_number = GetRandomNumber(hosts.GetN()-1);
        Ptr<Node> receiver = hosts.Get(reciever_number);

        // Random clients.
        std::vector<int> randomNumbers = GetRandomNumbers(num_of_flows, hosts.GetN()-1, reciever_number);  
        
        for(int i=0 ;i<randomNumbers.size() ;i++)
        {
            Ptr<Node> sender = hosts.Get(randomNumbers[i]);
            uint16_t port = basePort + i;

            // Install PacketSink on the receiver
            Address sinkAddress(InetSocketAddress(receiver->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal(), port));
            PacketSinkHelper sinkHelper("ns3::TcpSocketFactory", sinkAddress);
            ApplicationContainer sinkApp = sinkHelper.Install(receiver);
            sinkApp.Start(Seconds(1.0));
            sinkApp.Stop(Seconds(simulation_stop_time));
            sinkApps.push_back(sinkApp);

            // Install on off on the sender
            OnOffHelper onOffHelper("ns3::TcpSocketFactory", sinkAddress);
            onOffHelper.SetAttribute("DataRate", StringValue(app_data_rate)); // Sending rate
            onOffHelper.SetAttribute("PacketSize", UintegerValue(app_packet_size)); // Packet size (bytes)
            ApplicationContainer sourceApp = onOffHelper.Install(sender);
            sourceApp.Start(Seconds(2.0));
            sourceApp.Stop(Seconds(simulation_stop_time-1));


            if (i == 1) { // Track the first window
                Ptr<Application> app = sourceApp.Get(0); // Get the OnOffApplication instance
                Ptr<OnOffApplication> onOffApp = DynamicCast<OnOffApplication>(app);
                
                if (onOffApp) {
                    Simulator::Schedule(Seconds(2.01), [onOffApp]() {
                        Ptr<Socket> onOffSocket = onOffApp->GetSocket();
                        Ptr<TcpSocketBase> tcpSocketBase = DynamicCast<TcpSocketBase>(onOffSocket);
                        
                        if (tcpSocketBase) {
                            std::cout << "[DEBUG] Attaching cwnd tracer to connection 0 (Sender Node " << onOffApp->GetNode()->GetId() << ")" << std::endl;
                            Track_Flow(tcpSocketBase);  
                        } else {
                            std::cerr << "[ERROR] Could not cast OnOff socket to TcpSocketBase!" << std::endl;
                        }
                    });
                } else {
                    std::cerr << "[ERROR] Could not retrieve OnOffApplication instance!" << std::endl;
                }
            }
        }
    }
    
    else if(mode == 2)
    {
        // Random pairs of hosts to communicate.
        std::vector<Pair> random_pairs = GetRandomPairs(num_of_flows, hosts.GetN()-1);
        
        for(int i=0 ; i<random_pairs.size(); i++)
        {
            uint16_t port = basePort + i;
            Ptr<Node> receiver = hosts.Get(random_pairs[i].first);
            Ptr<Node> sender = hosts.Get(random_pairs[i].second);

            // Install PacketSink on the receiver
            Address sinkAddress(InetSocketAddress(receiver->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal(), port));
            PacketSinkHelper sinkHelper("ns3::TcpSocketFactory", sinkAddress);
            ApplicationContainer sinkApp = sinkHelper.Install(receiver);
            sinkApp.Start(Seconds(1.0));
            sinkApp.Stop(Seconds(simulation_stop_time));
            sinkApps.push_back(sinkApp);

            //  Install on off on the sender
            OnOffHelper onOffHelper("ns3::TcpSocketFactory", sinkAddress);
            onOffHelper.SetAttribute("DataRate", StringValue(app_data_rate)); // Sending rate
            onOffHelper.SetAttribute("PacketSize", UintegerValue(app_packet_size)); // Packet size (bytes)
            ApplicationContainer sourceApp = onOffHelper.Install(sender);
            sourceApp.Start(Seconds(2.0));
            sourceApp.Stop(Seconds(simulation_stop_time -1));


            if (i == 0) { // Track the first window
                Ptr<Application> app = sourceApp.Get(0); // Get the OnOffApplication instance
                Ptr<OnOffApplication> onOffApp = DynamicCast<OnOffApplication>(app);
                
                if (onOffApp) {
                    Simulator::Schedule(Seconds(2.01), [onOffApp]() {
                        Ptr<Socket> onOffSocket = onOffApp->GetSocket();
                        Ptr<TcpSocketBase> tcpSocketBase = DynamicCast<TcpSocketBase>(onOffSocket);
                        
                        if (tcpSocketBase) {
                            std::cout << "[DEBUG] Attaching cwnd tracer to connection 0 (Sender Node " << onOffApp->GetNode()->GetId() << ")" << std::endl;
                            Track_Flow(tcpSocketBase);  
                        } else {
                            std::cerr << "[ERROR] Could not cast OnOff socket to TcpSocketBase!" << std::endl;
                        }
                    });
                } else {
                    std::cerr << "[ERROR] Could not retrieve OnOffApplication instance!" << std::endl;
                }
            }
        }
    }
//-----------------------here ends the implementation of the connections between servers--------------------

   
//--------------------------flow monitor---------------------------
    int j=0;
    float AvgThroughput = 0;
    Time Delay = Seconds(0);
    uint32_t ReceivedPackets = 0;
    FlowMonitorHelper flowmon;
    float packet_loss_rate = 0;

    Ptr<FlowMonitor> monitor = flowmon.InstallAll ();

    Simulator::Schedule(Seconds(flowmonitor_start_time), [monitor]() {
    monitor->StartRightNow();
    });


    Simulator::Schedule(Seconds(simulation_stop_time + 10), [monitor]() {
    monitor->StopRightNow();
    });


    Simulator::Stop(Seconds(simulation_stop_time + 10));

    
    Simulator::Run();

    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmon.GetClassifier ());
    FlowMonitor::FlowStatsContainer stats = monitor->GetFlowStats ();

    for (auto iter = stats.begin (); iter != stats.end (); ++iter) 
    {
        Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (iter->first); 

        packet_loss_rate += (iter->second.txPackets - iter->second.rxPackets) * 100.0 / static_cast<double>(iter->second.txPackets);
        Delay = Delay + (iter->second.delaySum);
        ReceivedPackets = ReceivedPackets + (iter->second.rxPackets);
        AvgThroughput = AvgThroughput + (iter->second.rxBytes * 8.0/(iter->second.timeLastRxPacket.GetSeconds()-iter->second.timeFirstTxPacket.GetSeconds())/1024);
        j = j + 1;
        if (j == num_of_flows)
        {
            break;
        }
        
    }

    packet_loss_rate = packet_loss_rate / j;
    AvgThroughput = AvgThroughput/j;

    std::cout << "--------Total Results of the simulation----------"<<std::endl;
    std::cout << "Protocol: " << protocolName << std::endl;

    std::cout << "Average Throughput =" << AvgThroughput << "Kbps"<< std::endl;
    std::cout << "End to End Delay =" << Delay.GetSeconds()/ReceivedPackets<< std::endl;
    std::cout << "Average packet loss rate = " << packet_loss_rate <<"% " <<std::endl;

    //------------here starts the writing into the output file implementation---------

    std::string filename;

    filename = "scratch/results_" + protocolName + ".txt";

    std::ofstream outputFile(filename, std::ios::app);
    if (!outputFile.is_open()) 
    {
        std::cout << "Error: Could not open file output.txt for appending" << std::endl;
        return;
    }

    outputFile << "flows = " << num_of_flows << "   ";
    outputFile << "Throughput = " <<AvgThroughput << " Kbps   ";
    outputFile<< "Delay = " << Delay.GetSeconds()/ReceivedPackets<<"   ";
    outputFile << "time = " << simulation_stop_time << "   ";
    outputFile << "packet loss = " << packet_loss_rate <<"% " << "   K = " << k <<"    mode = "<< mode << std::endl;
    outputFile.close();
    //--------------here ends the writing into the output file implementation---------

    Simulator::Destroy();
    cwndFile.close();

    return;
}

//

int main(int argc, char *argv[])
{
    int counter = 0;
    int protocolNumber = atoi(argv[1]);

    SetTcpCongestionControl(protocolNumber);

    traceQueue = true;

    rttFile.open("scratch/" + protocolName + "_rtt.txt", std::ios::out);
    queueFile.open("scratch/" + protocolName + "_queue_size.txt", std::ios::out);
    cwndFile.open("scratch/" + protocolName + "_cwnd.txt", std::ios::out);

    if (protocolName == "TcpDctcp") {
        ecnFile.open("scratch/" + protocolName + "_ecn.txt", std::ios::out);
    }

    // runGit(4, "1Gbps", "10us", 1, 100, 0, 15, "100Mbps", 1000);
    // runGit(4, "1Gbps", "10us", 1, 120, 0, 15, "250Mbps", 1000);          // EXCEL
    // runGit(4, "1Gbps", "10us", 1, 120, 0, 30, "300Mbps", 1000);
    // OLD  runGit(4, "1Gbps", "10us", 1, 100, 0, 15, "100Mbps", 1000);
    // runGit(4, "1Gbps", "10us", 1, 120, 0, 30, "300Mbps", 1000);

    // 18.06.2025
    // runGit(4, "100Mbps", "10us", 1, 200, 0, 30, "500Mbps", 1000);   // הראשונים קורסים, בבר מציג 7.5%
    /////////////// runGit(4, "500Mbps", "10us", 1, 200, 0, 30, "500Mbps", 1000); // GOOD!!! Now small packets.
    
    // runGit(4, "500Mbps", "10us", 1, 20, 0, 30, "500Mbps", 70);  // Seems good! Still need to check
    // runGit(4, "500Mbps", "10us", 1, 200, 0, 30, "500Mbps", 70); // מלא זמן
    // runGit(4, "1Gbps", "10us", 1, 200, 0, 15, "500Mbps", 700);   // אין גודש
    
                                        /* 1 */
    // WOW!!!!                     0
    // runGit(4, "50Mbps", "10us", 1, 200, 0, 30, "500Mbps", 800); // GOOD!!!!!!!!
    // runGit(4, "50Mbps", "10us", 1, 200, 0, 30, "500Mbps", 1500); // WOW
    // runGit(4, "50Mbps", "10us", 1, 200, 0, 30, "500Mbps", 70); // Lot of time but OK

    // runGit(4, "40Mbps", "10us", 0, 200, 0, 30, "400Mbps", 700);  // Also good

    // runGit(4, "50Mbps", "10us", 1, 200, 0, 30, "500Mbps", 800); // GOOD!!!!!!!!


                                        /* 2 */
    // runGit(4, "1Mbps", "10us", 2, 400, 0, 30, "1Gbps", 700); // GOOD!!!
    // runGit(4, "1Mbps", "10us", 2, 400, 0, 30, "1Gbps", 1500); // Huge
    runGit(4, "1Mbps", "10us", 2, 400, 0, 30, "1Gbps", 70); // a lot of time. check if finish.

    // Small    61-70
    // Big      1400-1500
    // Average 800


    rttFile.close();
    queueFile.close();
    cwndFile.close();

    if (protocolName == "TcpDctcp") {
        ecnFile.close();
    }

    std::cout << std::endl;

    return 0;
}
