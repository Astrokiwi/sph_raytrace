#include <vector>
#include <memory>


struct SourceP_Type {
    double Pos[3];
    double brightness;
    double transmitting_area;
    double Hsml;
    double Mass;

};

//class absorbtree: public std::enable_shared_from_this<absorbtree> {
class absorbtree {
    public:
        absorbtree(double r_left[3], double size);
        //absorbtree(int level, int ixyz[3]);
        absorbtree(int level, int ixyz[3], double r_left[3], double size);
        
        void addPCoord(int ip,int ilevel,int ix,int iy,int iz);
        void head_addP(int ip, double r_p[3], double h);
        
        
        //SET_TYPE<int >* beam(double r1[3], double r2[3]);
        //SET_TYPE<int >* beam_level(int ilevel, double r1[3], double r2[3]);
        //void beam_level(int ilevel, double r1[3], double r2[3],std::vector<std::vector<int>*> *list_of_lists);
        //std::vector<std::vector<int>*>* beam(double r1[3], double r2[3]);

        void beam_level(int ilevel, double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize);
        void beam(double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize, int ThisTask);

        void beam_tree(double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize, int ThisTask);
        bool beam_subcells(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, double rtarget[3], unsigned char directionMap, std::vector<int>** lol, unsigned int* lolsize, int ThisTask);
        
        std::shared_ptr<absorbtree> getCellAt(int ilevel,int ix,int iy,int iz);
        std::shared_ptr<absorbtree> getCellAt(int ilevel,int ixyz[3]);
        std::vector<int>* getPList();
        
        void dump();
        void matrix_dump(int prefix);
        void sort();
        
        void murder();
        bool check(int N);
        
    private:
        inline unsigned char first_cell_hit(double tx0,double ty0,double tz0,double txm,double tym,double tzm);
        inline unsigned char choose_min(double x0,unsigned char i0, double x1, unsigned char i1, double x2, unsigned char i2);

        bool p_in_cell(double r[3]);
//        bool p_in_subcell(double r[3], int ipt[3]);

        int itreer(double r, int idim, int ilevel);
        double rtreei(int ii, int idim, int ilevel);

        std::vector<int> myP;

        double r_left[3];
        int ixyz[3];
        double size;
        int level;
        int maxlevel;

        std::shared_ptr<absorbtree> nodes[8];
        
//         std::ofstream dumpf; // FOR DEBUGGERING ONLY - REMOVE
//         bool mark; // FOR DEBUGGERING ONLY - REMOVE
        //PLIST_TYPE<int> *myP;
};

class maxlumtree {
    public:
        maxlumtree(double r_cent[3], double size);
        ~maxlumtree();
        void addp(double r_p[3], int ip, double temp);
        
        void dump(int ilevel);
        void dump();


//         std::vector<int>* getThreshList(double* allHotSourceP, double *r_p, int Ntot);
//         void propagateThreshList(double* allHotSourceP, double* r_p,std::vector<int>* threshList, int Ntot);
        void propagateThreshList(struct SourceP_Type* allHotSourceP, double* r_p, double receiving_area, int Ntot,int *threshList,unsigned int *nThresh);
        void getThreshList(struct SourceP_Type* allHotSourceP, double *r_p, double receiving_area, int Ntot, int* threshList, unsigned int* nThresh);

        double calcLuminosities(struct SourceP_Type* allHotSourceP);
        void buildNodeList();
        
        std::vector<maxlumtree*> *nodeList;
        double halfsize,size,size2;

        int listID;
        
        double r_cent[3];
        
        double luminosity; // sum of luminosity of all particles within

        int get_maxdepth();

    private:
        void placeInChildNode(double r_p[3], int ip, double temp);
        void propogateBuildNodeList(std::vector<maxlumtree*> *nodeList);

        int get_maxdepth(int maxlevel,int ilevel);
        
//         mintemptree **nodes;
        std::shared_ptr<maxlumtree> nodes[8];
        
        double r_p[3];
        int ip;
        double maxLum; // max luminosity of all particles within
        
        
};



class AGN_Kernel {
    public:
        AGN_Kernel();
        
        void calc_depth_optical_table(double depth, double mass, double loghmin, double loghmax, int nh);
        void set_L_flat(int L_in, int zL_in);


        double w(double x);
        double w(double r, double h);

        //double flat_w(double r, double h);
        double flat_w2(double r2, double h2);
        double half_flat_w2_truncated(double r2, double h2, double z);
        double skin_rad(double h, double opac);
    
    private:
        int L_flat,zL;
        double loghmin, loghmax;
        int nh;
        //std::vector<double> flattened_table;


        std::vector<double> skin_table;
        std::vector<double> flattened_table2;
        std::vector<std::vector<double>> fat_table2;
        
        void integrate_flat();
        
};


// separate out dust-dust heating into another file
void do_dust_dust_heating(bool alreadyDone[], int localOffset, int localNActive, int localIndex[], int nActiveTot,  double *tdiffwait,double *tdifftree,double *tdiffray);

double one_IR_tree_ray(std::shared_ptr<absorbtree> tree,
                       struct SourceP_Type allHotSourceP[],
                       bool alreadyDone[],double r_target[], double r_source[] /*r_source really*/,
                       double depthCrossSection,double d12_norm,double d12_norm2, double unextincted_flux);


// Table prototypes
// void setup_tables();
void setup_tables(
                    char AGNTabLabelFile[100], char AGNTabFile[100]
                    ,char AGNDustlessTabLabelFile[100], char AGNDustlessTabFile[100]
                    ,char AGNHighdenseTabLabelFile[100], char AGNHighdenseTabFile[100]
                  );
void agn_calc_heat_cool_p_raytraced_tab(int i);

#ifdef SOTON_AGN_PRES
double table_heat(int i, double *values);
#else
void table_heat(int i, double *values);
#endif  

// Raytracing prototypes
// std::shared_ptr<absorbtree> build_tree();
// void agn_optical_depths(double* r_agn, double *depths, std::shared_ptr<absorbtree> tree, bool agn_at_origin);
//void setup_raytracing();
//static inline double intersect_d2_nonorm(double *r1,double *r2,double *r3);

// global objects
// extern AGN_Kernel *agn_kernel;
// extern std::vector<int>** lol;
// extern const unsigned int maxlolsize;
