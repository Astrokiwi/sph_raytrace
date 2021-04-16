#include "rt_prototypes.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string> 
#include <limits>
#include <memory>

// #include "../allvars.h" // ONLY FOR DEBUGGERING - REMOVE!

/* Oct-tree for finding collisions between a ray and SPH particles.
Particles added to the tree are propagated down the oct-tree nodes
(creating new nodes along the way if possible) until they reach the 
shallowest level of refinement where cells are larger than the diameter
of the particle, i.e. 2*h. The particle is added to the particle-lists
for all cells at that refinement level that it intersects with, which
is up to eight cells.

This means that particles can be placed in branch nodes as well as in
leaf nodes. This stops particles with large softening lengths from being
placed in a very large number of cells - the max is 8.

A beam can be propagated through the tree. The beam is propagated through
all cells at each level of refinement in tern, returning a list of pointers
to particle-lists of cells intersected by the beam. This list will contain
duplicates.
*/

// std::ofstream dumpf;

// create head node
absorbtree::absorbtree(double r_left[3], double size) {
    for ( int ii=0 ; ii<3 ; ii++ ) {
        this->r_left[ii] = r_left[ii];
        this->ixyz[ii] = 0;
    }
    this->size = size;
    this->level = 0;
    
    this->maxlevel = 0;

}

// create leaf node
absorbtree::absorbtree(int level, int ixyz[3], double r_left[3], double size) {
    for ( int ii=0 ; ii<3 ; ii++ ) {
        this->r_left[ii] = r_left[ii];
        this->ixyz[ii] = ixyz[ii];
    }
    this->size = size;
    this->level = level;

}

bool absorbtree::check(int N) {
    if ( this->myP.size()>0 ) {
        for ( unsigned int ii=0 ; ii<this->myP.size() ; ii++ ) {
            if ( this->myP[ii]<0 || this->myP[ii]>=N ) {
                return false;
            }
        }
    }
    bool all_ok = true;
    for ( int ii=0 ; ii<8 ; ii++ ) {
        if ( nodes[ii] ) {
            if ( !nodes[ii]->check(N) ) {
                all_ok = false;
            }
        }
    }
    return all_ok;
}

void absorbtree::sort() {
    for ( int ii=0 ; ii<8 ; ii++ ) {
        if ( nodes[ii] ) {
            nodes[ii]->sort();
        }
    }
}
std::vector<int>* absorbtree::getPList() {
    return &(this->myP);
}

std::shared_ptr<absorbtree> absorbtree::getCellAt(int ilevel,int ix, int iy, int iz) {
    int ixyz[3];
    ixyz[0] = ix;
    ixyz[1] = iy;
    ixyz[2] = iz;
    return this->getCellAt(ilevel,ixyz);
}
    
std::shared_ptr<absorbtree> absorbtree::getCellAt(int ilevel,int ixyz[3]) {
    // determine which node
    int inode,dlevel;
    int ixyz_rel[3];
    dlevel = ilevel-this->level-1;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        ixyz_rel[ii] = (ixyz[ii]>>dlevel)-(this->ixyz[ii]<<1);
    }
    inode=ixyz_rel[0]+(ixyz_rel[1]<<1)+(ixyz_rel[2]<<2);

    if ( inode>=8 || inode<0 ) {
        return nullptr; // out of bounds
    }
    // check if node exists
    if ( !nodes[inode] ) {
        return nullptr; // no cell there
    }
    if ( ilevel==this->level+1 ) {
        return nodes[inode];
    } else {
        return nodes[inode]->getCellAt(ilevel,ixyz);
    }
}

void absorbtree::addPCoord(int ip,int ilevel,int ix,int iy,int iz) {
    if ( ilevel==this->level ) {
        myP.push_back(ip);
    } else {
        // determine which node
        int inode,dlevel;
        int ixyz[3],ixyz_rel[3];
        ixyz[0] = ix;
        ixyz[1] = iy;
        ixyz[2] = iz;
        dlevel = ilevel-this->level-1;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            ixyz_rel[ii] = (ixyz[ii]>>dlevel)-(this->ixyz[ii]<<1);
            if (ixyz_rel[ii]<0 || ixyz_rel[ii]>1) {
                if ( this->level==0 ) {
                    return; // just means cell is out of bounds
                } else {
                    std::cout << "CELL OUT OF BOUNDS AT REFINED LEVEL x<0 or x>1 " << ixyz_rel[ii] << std::endl;
                    std::cout << " " << ix << " " << iy << " " << iz << std::endl;
                    std::cout << " " << this->ixyz[0] << " " << this->ixyz[1] << " " << this->ixyz[2] << std::endl;
                    std::cout << ilevel << " " << this->level << " " <<dlevel << std::endl;
                    exit(0);
                }
            }

        }
        // ixyz_rel should be =0 or =1 for all three entries
        inode=ixyz_rel[0]+ixyz_rel[1]*2+ixyz_rel[2]*4;
        if ( inode>=8 || inode<0 ) {
            if ( this->level==0 ) {
                // This just means the cell is out of bounds.
                // This happens if the smoothing length overlaps the
                // edge of the tree of cells. If the tree contains all 
                // particles, then we can ignore this, because no ray
                // will pass outside this box when going between two
                // particles
                return;
            }
            // if ilevel!=0 then something is wrong
            exit(1);
            return;
        }


        // check if node exists
        if ( !nodes[inode] ) {
            int ixyz_node[3];
            double nodesize = size/2.;
            double r_left_node[3];
            for ( int ii=0 ; ii<3 ; ii++ ) {
                ixyz_node[ii] = (this->ixyz[ii]<<1)+ixyz_rel[ii];
                r_left_node[ii] = r_left[ii]+nodesize*ixyz_rel[ii];
            }
            nodes[inode] = std::make_shared<absorbtree>(this->level+1,ixyz_node,r_left_node,nodesize);
        }
        nodes[inode]->addPCoord(ip,ilevel,ix,iy,iz);
    }
}

// should only be called on the head node
// *requires* that particle is in grid
void absorbtree::head_addP(int ip, double r_p[3], double h) {
    int depth_cap = 8; // max is 256^2

    // goes in level where size of cell<hs
    int p_level;
    p_level = (int)floor(log2(size/h))-1; // +1 because we want 2h within each cell
    if ( p_level>depth_cap ) {
        p_level = depth_cap;
    }

    
    this->maxlevel = std::max(p_level,this->maxlevel);
    
    if ( p_level<=0 ) {
//         if (P[ip].ID==154731) { std::cout << "added to top level" << std::endl;}
        // add to myself, the top level
        myP.push_back(ip);
    } else {
        // add to lower levels
        // find which cells and level this particle should go
        double lsize = this->size/(1<<p_level); // cell size at target level
        int cent_coords[3]; // integer coordinates of left corner of cell at target level containing particle centre
        double diff_coords[3]; // displacement from corner of centre cell to particle centre
        int bump_over[3]; // which direction (if any) the particle overlaps with adjacent cels

        for ( int ii=0 ; ii<3 ; ii++ ) {
            cent_coords[ii] = (int)floor((r_p[ii]-r_left[ii])/lsize);
        }
        // add to centre cell
        this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]);
        
        // add to adjacent cells if particle overlaps. Loop over all 8 adjacent cells
        for ( int ii=0 ; ii<3 ; ii++ ) {
            diff_coords[ii] = r_p[ii]-r_left[ii]-lsize*cent_coords[ii];
            if ( h+diff_coords[ii]>=lsize ) {
                bump_over[ii] = 1;
            } else if ( diff_coords[ii]-h<0. ) {
                bump_over[ii] = -1;
            } else {
                bump_over[ii] = 0;
            }
        }
        // place in all adjacent cells
        if ( bump_over[0]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1],cent_coords[2]);
            if ( bump_over[1]!=0 ) {
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]);
                this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1]+bump_over[1],cent_coords[2]);
                if ( bump_over[2]!=0 ) {
                    this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
                    this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                }
            } else if ( bump_over[2]!=0 ) {
                this->addPCoord(ip,p_level,cent_coords[0]+bump_over[0],cent_coords[1],cent_coords[2]+bump_over[2]);
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
            }
        } else if ( bump_over[1]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]);
            if ( bump_over[2]!=0 ) {
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1]+bump_over[1],cent_coords[2]+bump_over[2]);
                this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
            }
        } else if ( bump_over[2]!=0 ) {
            this->addPCoord(ip,p_level,cent_coords[0],cent_coords[1],cent_coords[2]+bump_over[2]);
        }
    }
}

// return the set of all particles in cells that intersect the beam
// only called on head
void absorbtree::beam(double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize,int ThisTask) {

   if ( this->level!=0 ) return;

// check if beam passes through *any* cell at all
    bool passes_through = true;
    for ( int ii=0 ; ii<3 ; ii++ ) {
        if ( r1[ii]<r_left[ii] && r2[ii]<r_left[ii] ) passes_through = false;
        if ( r1[ii]>r_left[ii]+size && r2[ii]>r_left[ii]+size ) passes_through = false;
    }
    if ( !passes_through ) return;

  beam_tree(r1,r2,lol,lolsize,ThisTask);
}

void absorbtree::beam_tree(double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize, int ThisTask) {
    double ddr[3];

    double t0[3],t1[3];
    
    unsigned char directionMap = 0;
    double r_start[3];
    
    // calculate initial directions
    for ( int ii=0 ; ii<3 ; ii++ ) {
        ddr[ii] = r2[ii]-r1[ii];
        if ( ddr[ii]==0. ) {
            ddr[ii]=size*1.e-9;
        }


        if (ddr[ii]<0.) {
            r_start[ii] = this->size-r1[ii]+2.*this->r_left[ii];
            ddr[ii]=-ddr[ii];
            directionMap |= 1<<ii;
        } else {
            r_start[ii] = r1[ii];
        }
    
        t0[ii] = (this->r_left[ii]-r_start[ii])/ddr[ii];
        t1[ii] = (this->r_left[ii]+this->size-r_start[ii])/ddr[ii];
    }

    
    beam_subcells(t0[0],t0[1],t0[2],t1[0],t1[1],t1[2],r2,directionMap,lol,lolsize,ThisTask);
}

bool absorbtree::p_in_cell(double r[3]) {
    for ( int ii = 0 ; ii<3 ; ii++ ) {
        if ( r[ii]<this->r_left[ii] || r[ii]>=this->r_left[ii]+this->size ) {
            return false;
        }
    }
    return true;
}

inline unsigned char absorbtree::first_cell_hit(double tx0,double ty0,double tz0,double txm,double tym,double tzm) {
    unsigned int plane;
    if ( tx0>ty0 ) {
        if ( tz0>tx0 )
            plane = 2;
        else
            plane = 0;
    } else {
        if ( tz0>ty0 )
            plane = 2;
        else
            plane = 1;
    }
    
    unsigned char node_out = 0;
    switch(plane) {
        case 0:
            if ( tym<tx0 ) node_out|=2;
            if ( tzm<tx0 ) node_out|=4;
            break;
        case 1:
            if ( txm<ty0 ) node_out|=1;
            if ( tzm<ty0 ) node_out|=4;
            break;
        case 2:
            if ( txm<tz0 ) node_out|=1;
            if ( tym<tz0 ) node_out|=2;
            break;
    }
    
    return node_out;
}

inline unsigned char absorbtree::choose_min(double x0,unsigned char i0, double x1, unsigned char i1, double x2, unsigned char i2) {
    if ( x0<x1 ) {
        if ( x2<x0 ) return i2;
        else return i0; 
    } else {
        if ( x2<x1 ) return i2;
        else return i1;
    }
}

// Revelles, Ureña, Lastra 2000 algorithm, but with different indexing for nodes (xyz instead of zyx)
bool absorbtree::beam_subcells(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1, double rtarget[3], unsigned char directionMap, std::vector<int>** lol, unsigned int* lolsize, int ThisTask) {
    if ( this->myP.size()>0 ) {
        lol[*lolsize] = &this->myP;
        (*lolsize)++;
    }

    if ( tx1<0. ) return false;//p_in_cell(rtarget);
    if ( ty1<0. ) return false;//p_in_cell(rtarget);
    if ( tz1<0. ) return false;//p_in_cell(rtarget);

    double txm,tym,tzm;
    unsigned char inode;
    bool hitp=false;
    unsigned char jnode;

    txm = .5*(tx0+tx1);
    tym = .5*(ty0+ty1);
    tzm = .5*(tz0+tz1);
    
    inode = first_cell_hit(tx0,ty0,tz0,txm,tym,tzm);
///////////////
    do {
        jnode = inode^directionMap;
//         std::cout << (int)inode << " " << (int)jnode << std::endl;
        switch(inode) {
            case 0:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(tx0,ty0,tz0,txm,tym,tzm,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(txm,1,tym,2,tzm,4);
                break;
            case 1:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(txm,ty0,tz0,tx1,tym,tzm,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(tx1,8,tym,3,tzm,5);
                break;
            case 2:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(tx0,tym,tz0,txm,ty1,tzm,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(txm,3,ty1,8,tzm,6);
                break;
            case 3:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(txm,tym,tz0,tx1,ty1,tzm,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(tx1,8,ty1,8,tzm,7);
                break;
            case 4:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(tx0,ty0,tzm,txm,tym,tz1,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(txm,5,tym,6,tz1,8);
                break;
            case 5:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(txm,ty0,tzm,tx1,tym,tz1,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(tx1,8,tym,7,tz1,8);
                break;
            case 6:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(tx0,tym,tzm,txm,ty1,tz1,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = choose_min(txm,7,ty1,8,tz1,8);
                break;
            case 7:
                if ( nodes[jnode] ) hitp=nodes[jnode]->beam_subcells(txm,tym,tzm,tx1,ty1,tz1,rtarget,directionMap,lol,lolsize,ThisTask);
                inode = 8;
                break;
        }
    } while (inode<8 && !hitp);

    if ( hitp ) {
        return true;
    } else {
        return p_in_cell(rtarget);
    }
}

// find grid coordinate at this level
int absorbtree::itreer(double r, int idim, int ilevel) {
    return (int)floor(((r-this->r_left[idim])/this->size)*(1<<ilevel));
}

// find "left" corner of cell at given coordinate
double absorbtree::rtreei(int ii, int idim, int ilevel) {
    return (ii*this->size)/(1<<ilevel)+this->r_left[idim];
}

// return the set of all particles in cells at this particular level of refinement that intersect the beam
// only called on head
void absorbtree::beam_level(int ilevel, double r1[3], double r2[3],std::vector<int>** lol, unsigned int* lolsize) {
   
   if ( this->level!=0 ) return;

   int idest[3];
   for ( int ii=0 ; ii<3 ; ii++ ) {
       idest[ii] = itreer(r2[ii],ii,ilevel);
    }

    
    // Move through cells at this level
    double ddr[3],liner[3],dt[3];
    int iline[3],idir[3];
    int dir_choice;
    bool ingrid=true,hitp=false;
    
    int L = 1<<ilevel;
    
    // calculate initial directions
    for ( int ii=0 ; ii<3 ; ii++ ) {
        ddr[ii] = r2[ii]-r1[ii];
        liner[ii] = r1[ii];
        iline[ii] = itreer(r1[ii],ii,ilevel);
        if ( ddr[ii]>0. ) {
            idir[ii] = 1;
        } else if ( ddr[ii]<0. ) {
            idir[ii] = -1;
        } else {
            idir[ii] = 0;
        }
    }
    
    if ( ilevel==0 ) {
        lol[*lolsize] = this->getPList();
        (*lolsize)++;
    } else {
        std::shared_ptr<absorbtree> cell0 = this->getCellAt(ilevel,iline[0],iline[1],iline[2]);
        if ( cell0 ) {
            lol[*lolsize] = cell0->getPList();
            (*lolsize)++;
        }
    }
    
    // Loop through line
    // This version is inefficient and loops through lots of empty cells, especially at greatest refinement level!
    // TODO: make less dumb
    while(ingrid && !hitp) {

        // calculate "dt" to each edge
        // go in the direction with the lowest "dt"
        // this will currently fail if ddr[ii]=0. for any direction, so we'll have to fix that
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( idir[ii]==1 )         dt[ii] = (rtreei(iline[ii]+1,ii,ilevel)-liner[ii])/ddr[ii]; // "right" side of this cell is "left" cell of next cell
            else if ( idir[ii]==-1 )   dt[ii] = (rtreei(iline[ii],ii,ilevel)-liner[ii])/ddr[ii]; // "left" side of this cell
            else dt[ii] = std::numeric_limits<double>::max();//(rtreei(iline[ii]+idir[ii],ii,ilevel)-liner[ii])/ddr[ii];
        }
        dir_choice = 0;
        for ( int ii=1 ; ii<3 ; ii++ ) {
            if ( dt[ii]<dt[dir_choice] ) {
                dir_choice = ii;
            }
        }
        
        // jump to the next cell
        iline[dir_choice]+=idir[dir_choice];
        
        // update line location
        for ( int ii=0 ; ii<3 ; ii++ ) {
            liner[ii]+=dt[dir_choice]*ddr[ii];
        }
        
        // check for finishing conditions
        ingrid = true;
        hitp = true;
        for ( int ii=0 ; ii<3 ; ii++ ) {
            if ( iline[ii]<0 || iline[ii]>=L ) ingrid=false;
            if ( iline[ii]!=idest[ii] ) hitp = false;
        }
        
        // add the list of particles
        if ( ingrid ) {
            std::shared_ptr<absorbtree> hitCell = this->getCellAt(ilevel,iline[0],iline[1],iline[2]);
            if ( hitCell ) {
                lol[*lolsize] = hitCell->getPList();
                (*lolsize)++;
            }
        }
    }
}



void absorbtree::dump() {
    for ( int ii=0; ii<this->level; ii++ ) {
        std::cout << ".";
    }
    std::cout << myP.size() << ":";
    for ( auto const& pindex: myP ) {
        std::cout << pindex << " ";
    }
    std::cout << std::endl;
    for ( int ii=0; ii<8; ii++ ) {
        if ( nodes[ii] ) {
            nodes[ii]->dump();
        } else {
            for ( int ii=0; ii<this->level; ii++ ) {
                std::cout << ".";
            }
            std::cout << "-" << std::endl;
        }
    }
}

void absorbtree::matrix_dump(int prefix) {
    for ( int ilevel = 0 ; ilevel<=this->maxlevel ; ilevel++ ) {
        std::ofstream f;
        f.open("treedump"+std::to_string(prefix)+"_"+std::to_string(ilevel)+".dat");
        int L = 1<<ilevel;
        for ( int ix=0 ; ix<L ; ix++ ) {
            for ( int iy = 0 ; iy<L ; iy++ ) {
                int sizeTot = 0;
                for ( int iz=0 ; iz<L ; iz++ ) { 
                    std::shared_ptr<absorbtree> thisCell = this->getCellAt(ilevel,ix,iy,iz);
                    if ( !thisCell ) {
                        sizeTot+=thisCell->getPList()->size();
                    }
                }
                f << sizeTot << " ";
            }
            f << std::endl;
        }
        f.close();
    }
}






