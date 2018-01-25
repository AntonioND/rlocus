/*
    Root Locus 0.2, Copyright (C) 2011-2012 Antonio Niño Díaz

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <hpgcc49.h>
#include <hpctype.h>
#include <hpgraphics.h>
#include <hpmath.h>
#include <hpstdio.h>

//-----------------------------------------------------------

void error(char * txt)
{
    //hpg_set_mode_mono(0);
    hpg_clear();
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text(txt,0,0);
}

#define SCR_X 131
#define SCR_Y 80

#define MAX_DEGREE (12+1)

//-------------------------------------------------------------------------------------

#include "PRF/header.h"

/*dcomplex *p,     coefficient vector of the original polynomial
           *pred,  coefficient vector of the deflated polynomial
           *root;  determined roots
  int      *n;     the highest exponent of the original polynomial
  double   *maxerr; max. error of all determined roots
  unsigned char flag;  flag = TRUE  => complex coefficients
                     flag = FALSE => real    coefficients
*/
unsigned char null(dcomplex *p,dcomplex *pred,int *n,dcomplex *root, double *maxerr,unsigned char flag);

dcomplex poly[MAX_DEGREE], pred[MAX_DEGREE], roots[MAX_DEGREE];

int solve_poly(double * a, int n, double * sr,double * si)
{
    int max_deg = n;
    double max_err;

    int i;
    for(i = 0; i <= n; i++) { poly[i].i = 0.0; poly[i].r = a[n-i]; }

    int reduced = 0;

    while( (poly[0].r == 0) && (poly[0].i == 0) )
    {
        int i;
        for(i = 0; i < n; i++) { poly[i].r = poly[i+1].r; poly[i].i = poly[i+1].i; }

        max_deg --;
        reduced = 1;
    }

    int reduced_n = max_deg;

    if(null(poly,pred,&max_deg,roots,&max_err,TRUE))
    {
        if(reduced) //return "x=0"
            for(i = 0; i < n; i++) { sr[i] = 0; si[i] = 0; }
        else
            return 0; //error, no roots
    }

    for(i = 0; i < reduced_n; i++) { sr[i] = roots[i].r; si[i] = roots[i].i; }
    for( ; i < n; i++) { sr[i] = 0; si[i] = 0; }

    return n;
}

//-------------------------------------------------------------------------------------

double num_coefs[MAX_DEGREE]; //COEFS[0] == x^n [...] COEFS[DEGREE] == x^0 || ( x^n ...  x^2  x^1  x^0 )
int NUM_DEGREE; // DEGREE = n
double den_coefs[MAX_DEGREE]; //COEFS[0] == x^n [...] COEFS[DEGREE] == x^0 || ( x^n ...  x^2  x^1  x^0 )
int DEN_DEGREE; // DEGREE = n

int is_valid_num(char c)
{
    return ( isdigit(c) || (c == '.') || (c == ',') || (c == 'E') || (c == '-') );
}

//----------------------
double saved_num_coefs[MAX_DEGREE], saved_den_coefs[MAX_DEGREE];
int saved_NUM_DEGREE, saved_DEN_DEGREE;
int saved = 0;

void save_parsed_polys(void)
{
    saved = 1;
    saved_NUM_DEGREE = NUM_DEGREE; saved_DEN_DEGREE = DEN_DEGREE;
    int i;
    for(i = 0; i < MAX_DEGREE; i++)
        { saved_num_coefs[i] = num_coefs[i]; saved_den_coefs[i] = den_coefs[i]; }
}

void restore_parsed_polys(void)
{
    if(saved == 0) return;
    NUM_DEGREE = saved_NUM_DEGREE; DEN_DEGREE = saved_DEN_DEGREE;
    int i;
    for(i = 0; i < MAX_DEGREE; i++)
        { num_coefs[i] = saved_num_coefs[i]; den_coefs[i] = saved_den_coefs[i]; }
}

//----------------------

int parse_poly(void)
{
    unsigned rpl_stack_bias = sat_stack_init();
/*
    sat_stack_push_string("{ 1 }");
    sat_stack_push_string("{ 1 0 }");
*/
/*
    sat_stack_push_string("{ 1 2 3  }");
    sat_stack_push_string("{ 1 3 2 0 0 }");
*/
/*
    sat_stack_push_string("{ 1.0 1.0}");
    //sat_stack_push_string("{ 1.0 -2.0 1.0 -2.0 }");
    sat_stack_push_string("{ 1.0 5.0 6.0 }");
*/
/*
    sat_stack_push_string("{ 1 1}");
    sat_stack_push_string("{ 1 -1 0 }");
*/
/*
    sat_stack_push_string("{ 1 }");
    sat_stack_push_string("{ 1 -0.9 0.25 0.225 }");
*/
/*
    sat_stack_push_string("{ 1 0.3 }");
    sat_stack_push_string("{ 1 -0.5 0.5 }");
*/
/*
    sat_stack_push_string("{ 1.0 3.5}");
    sat_stack_push_string("{ 1.0 -2.0 1.0 -2.0 }");
*/
/*
    sat_stack_push_string("{ 1.0 3.5}");
    sat_stack_push_string("{ 1 17 94 168 }");
*/
/*
    sat_stack_push_string("{ 1.0 }");
    sat_stack_push_string("{ 1 9 82 192 }");
*/
/*
    sat_stack_push_string("{ 1 1}");
    sat_stack_push_string("{ 1 -3 1 }");
*/
/*
    sat_stack_push_string("{ 1 3}");
    sat_stack_push_string("{ 1 3 2 }");
*/

    if(sat_stack_depth()<2){ //check for enough items
        error("Error: Too Few Arguments!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    SAT_STACK_ELEMENT first;
    SAT_STACK_ELEMENT second;
    sat_get_stack_element(1,&first);
    sat_get_stack_element(1,&second);
    if( (first.prologue!=SAT_DOCSTR) || (second.prologue!=SAT_DOCSTR) ) {
        error("Error: Invalid Argument!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    // DENOMINATOR
    //-------------

    char str[100];
    sat_pop_str(str);

    int l = strlen(str);
    if( (str[0] != '{') || (str[l-1] != '}') )
    {
        error("Error: Invalid Argument!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    str[0] = ' ';
    str[l-1] = ' ';

    memset((void*)den_coefs,0,sizeof(den_coefs));
    DEN_DEGREE = 0;

    int i = 0;
    while(i < l) //get number of coefficients
    {
        char c = str[i++];
        if( is_valid_num(c) )
        {
            DEN_DEGREE++;
            while(1)
            {
                c = str[i++];
                if( !is_valid_num(c) ) break;
                if(i >= l) //can't end in number - not a list
                {
                    error("Error: Invalid Argument!");
                    sat_stack_exit(rpl_stack_bias);
                    return 1;
                }
            }
        }
        else if(c != ' ')
        {
            error("Error: Invalid Argument!");
            sat_stack_exit(rpl_stack_bias);
            return 1;
        }
    }

    if(DEN_DEGREE > MAX_DEGREE)
    {
        error("Error: Too big polynomial!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    i = 0;
    char * strparse = str;
    strparse ++;
    while(*strparse)
    {
        char c = *strparse;
        if(c != ' ') break;
        strparse++;
    }
    while(i < DEN_DEGREE) //DEGREE IS THE NUMBER OF COEFFICIENTS
    {
        if(sscanf(strparse,"%lf",&(den_coefs[i])) != 1)
        {
            error("Error: Invalid Argument!");
            sat_stack_exit(rpl_stack_bias);
            return 1;
        }

        i++;

        while(*strparse)
        {
            char c = *strparse;
            if(c != ' ') break;
            strparse++;
        }

        while(*strparse)
        {
            if(!is_valid_num(*strparse)) break;
            strparse++;
        }

        while(*strparse)
        {
            char c = *strparse;
            if(c != ' ') break;
            strparse++;
        }
    }

    DEN_DEGREE--; //DEGREE WAS THE NUMBER OF COEFFICIENTS

    // NUMERATOR
    //-----------

    sat_pop_str(str);

    l = strlen(str);
    if( (str[0] != '{') || (str[l-1] != '}') )
    {
        error("Error: Invalid Argument!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    str[0] = ' ';
    str[l-1] = ' ';

    memset((void*)num_coefs,0,sizeof(num_coefs));
    NUM_DEGREE = 0;

    i = 0;
    while(i < l) //get number of coefficients
    {
        char c = str[i++];
        if( is_valid_num(c) )
        {
            NUM_DEGREE++;
            while(1)
            {
                c = str[i++];
                if( !is_valid_num(c) ) break;
                if(i >= l) //can't end in number - not a list
                {
                    error("Error: Invalid Argument!");
                    sat_stack_exit(rpl_stack_bias);
                    return 1;
                }
            }
        }
        else if(c != ' ')
        {
            error("Error: Invalid Argument!");
            sat_stack_exit(rpl_stack_bias);
            return 1;
        }
    }

    if(NUM_DEGREE > MAX_DEGREE)
    {
        error("Error: Too big polynomial!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    i = 0;
    strparse = str;
    strparse ++;
    while(*strparse)
    {
        char c = *strparse;
        if(c != ' ') break;
        strparse++;
    }
    while(i < NUM_DEGREE) //DEGREE IS THE NUMBER OF COEFFICIENTS
    {
        if(sscanf(strparse,"%lf",&(num_coefs[i])) != 1)
        {
            error("Error: Invalid Argument!");
            sat_stack_exit(rpl_stack_bias);
            return 1;
        }

        i++;

        while(*strparse)
        {
            char c = *strparse;
            if(c != ' ') break;
            strparse++;
        }

        while(*strparse)
        {
            if(!is_valid_num(*strparse)) break;
            strparse++;
        }

        while(*strparse)
        {
            char c = *strparse;
            if(c != ' ') break;
            strparse++;
        }
    }

    NUM_DEGREE--; //DEGREE WAS THE NUMBER OF COEFFICIENTS

    if((NUM_DEGREE < 0) || (DEN_DEGREE < 0))
    {
        error("Error: Empty list!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    if( (num_coefs[0] == 0) || (den_coefs[0] == 0) )
    {
        error("Error: First coefficient is 0!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    if((NUM_DEGREE < 1) && (DEN_DEGREE < 1))
    {
        error("Error: Input is constant!");
        sat_stack_exit(rpl_stack_bias);
        return 1;
    }

    sat_stack_exit(rpl_stack_bias);

    save_parsed_polys();

    return 0;
}

double COEFS[MAX_DEGREE];
int DEGREE;
void multiply_poly_add(double k)
{
    if(NUM_DEGREE > DEN_DEGREE) DEGREE = NUM_DEGREE;
    else DEGREE = DEN_DEGREE;

    int i;
    for(i = 0; i < (DEGREE+1); i++)
    {
        double a = 0;

        if( NUM_DEGREE < DEGREE )
        {
            int diff = DEN_DEGREE-NUM_DEGREE;

            if(i >= diff)
                a += num_coefs[i-diff]*k;

            a += den_coefs[i];
        }
        else
        {
            int diff = NUM_DEGREE-DEN_DEGREE;

            if(i >= diff)
                a += den_coefs[i-diff];

            a += num_coefs[i]*k;
        }

        COEFS[i] = a;
    }
}

//use num_coefs/NUM_DEGREE || den_coefs/DEN_DEGREE
void multiply_poly_regulator(double * coefs, int * degree, double a) // multiply coefs by (x+a)
{
    double new_coefs[MAX_DEGREE];

    int i;
    for(i = 0; i <= (1+*degree); i++)
    {
        if(i == 0) new_coefs[0] = coefs[0];
        else if(i <= *degree) new_coefs[i] = coefs[i] + coefs[i-1]*a;
        else new_coefs[i] = coefs[i-1]*a;
    }

    *degree=1+*degree;

    for(i = 0; i < (1+*degree); i++) coefs[i] = new_coefs[i];
}

void multiply_poly_const(double * coefs, int degree, double k) // multiply coefs by k
{
    int i;
    for(i = 0; i < (1+degree); i++) coefs[i] *= k;
}

//--------------------------------------------------------------

int get_double(char * show_string, double * returned_value) //returns 1 if OK
{
    char string[30+1];
    int count = 0;

    hpg_clear();
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text(show_string,0,0);

    while(1)
    {
        if(keyb_isAnyKeyPressed())
        {
            if(keyb_isON()) return 0;

            if(keyb_isKeyPressed(0,6)) { while(keyb_isKeyPressed(0,6)) sys_LCDSynch(); break; } //ENTER

            if(count < 30)
            {
                struct { int c, r; char value; } info[15] = {
                    {3,6,'0'}, {3,5,'1'}, {2,5,'2'}, {1,5,'3'}, {3,4,'4'}, {2,4,'5'}, {1,4,'6'}, {3,3,'7'},
                    {2,3,'8'}, {1,3,'9'}, {2,6,'.'}, {4,2,'E'}, {5,4,'E'}, {0,4,'-'}, {3,2,'-'}
                };

                int i = 0;
                while(i < 15)
                {
                    if(keyb_isKeyPressed(info[i].c,info[i].r))
                    {
                        string[count++] = info[i].value;
                        while(keyb_isKeyPressed(info[i].c,info[i].r)) sys_LCDSynch();
                        break;
                    }
                    i++;
                }
            }

            if(count > 0) if(keyb_isKeyPressed(0,0)) { count--; while(keyb_isKeyPressed(0,0)) sys_LCDSynch(); }

            string[count] = '\0';

            hpg_set_color(hpg_stdscreen,HPG_COLOR_WHITE);
            hpg_fill_rect(0,8,130,16);
            hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
            hpg_draw_text(string,0,8);
        }

        sys_LCDSynch();
    }

    if(count == 0) return 0;

    double res;
    int ok = sscanf(string,"%lf",&res);

    *returned_value = res;

    return ok;
}

//--------------------------------------------------------------

int SYSTEM_IS_CONTINOUS = 1;

double get_x(double i) //i MUSTN'T BE 100.00 ( draw_paths uses i < 100.00 )
{
    if(i < 50.0) return pow(i/50.0,2);
    return pow(1.0/(2.0-(i/50.0)),2);
}

double scale_factor;
double axes_scale;
int r_offset;

void plot(double xx, double yy)
{
    int x = xx*scale_factor+(double)r_offset;
    int y = -yy*scale_factor+((double)SCR_Y/2.0);
    if(x < 0 || x >= SCR_X) return;
    if(y < 0 || y >= SCR_Y) return;
    hpg_draw_pixel(x,y);
}

int test_plot(double xx, double yy)
{
    int x = xx*scale_factor+(double)r_offset;
    int y = -yy*scale_factor+((double)SCR_Y/2.0);
    if(x < 0 || x >= SCR_X) return 0;
    if(y < 0 || y >= SCR_Y) return 0;
    return 1;
}

void plot_add(double xx, double yy, int xadd, int yadd)
{
    int x = xx*scale_factor+(double)r_offset;
    int y = -yy*scale_factor+((double)SCR_Y/2.0);
    x+=xadd;
    y+=yadd;
    if(x < 0 || x >= SCR_X) return;
    if(y < 0 || y >= SCR_Y) return;
    hpg_draw_pixel(x,y);
}

void auto_scale(void)
{
    double maxr = 0.1, minr = -0.1, maxi = 0.1, mini = -0.1;
    double defaultval = 0.1;

    while(1)
    {
        int modified = 0;

        double wr_n[NUM_DEGREE],wi_n[NUM_DEGREE];
        int numr_n = solve_poly(num_coefs,NUM_DEGREE,wr_n,wi_n);
        int i; for(i = 0; i < numr_n; i++)
        {
            if(wr_n[i] > maxr) { maxr = wr_n[i]; modified = 1; }
            if(wr_n[i] < minr) { minr = wr_n[i]; modified = 1; }
            if(wi_n[i] > maxi) { maxi = wi_n[i]; modified = 1; }
            if(wi_n[i] < mini) { mini = wi_n[i]; modified = 1; }
        }

        double wr_d[DEN_DEGREE],wi_d[DEN_DEGREE];
        int numr_d = solve_poly(den_coefs,DEN_DEGREE,wr_d,wi_d);
        for(i = 0; i < numr_d; i++)
        {
            if(wr_d[i] > maxr) { maxr = wr_d[i]; modified = 1; }
            if(wr_d[i] < minr) { minr = wr_d[i]; modified = 1; }
            if(wi_d[i] > maxi) { maxi = wi_d[i]; modified = 1; }
            if(wi_d[i] < mini) { mini = wi_d[i]; modified = 1; }
        }

        if(modified) break;
        else
        {
            maxr /= 10.0;
            minr /= 10.0;
            maxi /= 10.0;
            mini /= 10.0;
            defaultval /= 10.0;

            if(defaultval == 0.0)
            {
                //everything is in x=0? ¬_¬
                //PD: this is here because of Susana V.
                maxr = 1;
                minr = -1;
                maxi = 1;
                mini = -1;
                defaultval = 1;
                break;
            }
        }
    }

    if(maxr == defaultval)
    {
        if(minr != defaultval)
        {
            maxr = fabs(minr)/4.0;
        }
    }
    if(minr == defaultval)
    {
        if(maxr != defaultval)
        {
            minr = -fabs(maxr)/4.0;
        }
    }

    maxr *= 2.0;
    minr *= 2.0;
    maxi *= 2.0;
    mini *= 2.0;

    double r_range = maxr-minr;
    double i_range = maxi-mini;

    r_offset = (-minr*(double)SCR_X)/r_range;

    if(r_range/i_range > ((double)SCR_X/(double)SCR_Y)) i_range = (r_range * (double)SCR_Y) / (double)SCR_X;
    else r_range = (i_range * (double)SCR_X) / (double)SCR_Y;

    scale_factor = ((double)SCR_X)/r_range;
}

void draw_axes(void)
{
    hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_5);
    hpg_draw_line(0,SCR_Y/2,SCR_X-1,SCR_Y/2);
    hpg_draw_line(r_offset,0,r_offset,SCR_Y);

    double increment = pow(10,floor(log10(((double)SCR_X)/scale_factor)));
    while(1) //check Y axis
    {
        if(test_plot(0,increment)) increment *= 10.0;
        else break;
    }
    increment /= 10.0;

    axes_scale = increment;

    if(SYSTEM_IS_CONTINOUS == 0)
    {
        //Draw circle with radius = 1.0
        double a, inc_circ = (M_PI*2.0)/(scale_factor*100.0);
        for(a = 0.0; a < 2.0*M_PI; a += inc_circ)
            plot(cos(a),sin(a));
    }

    //X AXIS
    //------

    double xx = 0;
    while(1)
    {
        int x = xx*scale_factor+(double)r_offset;
        if(x < 0 || x >= SCR_X) break;
        hpg_draw_pixel(x,(SCR_Y/2) -1);
        xx -= increment;
    }

    xx = 0;
    while(1)
    {
        int x = xx*scale_factor+(double)r_offset;
        if(x < 0 || x >= SCR_X) break;
        hpg_draw_pixel(x,(SCR_Y/2) -1);
        xx += increment;
    }

    //Y AXIS
    //------

    int add_y_axis;
    if((r_offset+1) < SCR_X) add_y_axis = 1;
    else add_y_axis = -1;

    double yy = 0;
    while(1)
    {
        int x = r_offset+add_y_axis;
        int y = yy*scale_factor+((double)SCR_Y/2.0);
        if(y < 0 || y >= SCR_Y) break;
        hpg_draw_pixel(x,y);
        yy -= increment;
    }

    yy = 0;
    while(1)
    {
        int x = r_offset+add_y_axis;
        int y = yy*scale_factor+((double)SCR_Y/2.0);
        if(y < 0 || y >= SCR_Y) break;
        hpg_draw_pixel(x,y);
        yy += increment;
    }
}

void draw_singularities(void)
{
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);

    //ZEROS
    double wr_n[NUM_DEGREE],wi_n[NUM_DEGREE];
    int numr_n = solve_poly(num_coefs,NUM_DEGREE,wr_n,wi_n);
    int i; for(i = 0; i < numr_n; i++)
    {
        int x = wr_n[i]*scale_factor+(double)r_offset;
        int y = wi_n[i]*scale_factor+((double)SCR_Y/2.0);
        //no need to test, auto scale should handle it
        hpg_draw_pixel(x+2,y);
        hpg_draw_pixel(x+1,y+1);
        hpg_draw_pixel(x,y+2);
        hpg_draw_pixel(x-1,y+1);
        hpg_draw_pixel(x-2,y);
        hpg_draw_pixel(x-1,y-1);
        hpg_draw_pixel(x,y-2);
        hpg_draw_pixel(x+1,y-1);
    }

    //POLES
    double wr_d[DEN_DEGREE],wi_d[DEN_DEGREE];
    int numr_d = solve_poly(den_coefs,DEN_DEGREE,wr_d,wi_d);
    for(i = 0; i < numr_d; i++)
    {
        int x = wr_d[i]*scale_factor+(double)r_offset;
        int y = wi_d[i]*scale_factor+((double)SCR_Y/2.0);
        hpg_draw_pixel(x+2,y+2);
        hpg_draw_pixel(x+1,y+1);
        hpg_draw_pixel(x,y);
        hpg_draw_pixel(x-1,y-1);
        hpg_draw_pixel(x-2,y-2);
        hpg_draw_pixel(x-1,y+1);
        hpg_draw_pixel(x-2,y+2);
        hpg_draw_pixel(x+1,y-1);
        hpg_draw_pixel(x+2,y-2);
    }
}

void draw_paths(void)
{
    int speedup = 0;

    int lastpercent = 0;

    double i = 1.0;
    while(i < 100.00)
    {
        multiply_poly_add(get_x(i));
        hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_10);
        double wr[DEGREE],wi[DEGREE];
        int numr = solve_poly(COEFS,DEGREE,wr,wi);

        if(numr > 0)
        {
            int j;
            for(j=0; j<numr; j++) plot(wr[j],wi[j]);
        }

        if(lastpercent != (int)i)
        {
            lastpercent = (int)i;
            char str[20];
            sprintf(str,"%d%%",lastpercent);
            hpg_set_color(hpg_stdscreen,HPG_COLOR_WHITE);
            hpg_fill_rect(0,0,20,5);
            hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
            hpg_draw_text(str,0,0);
        }
        if(keyb_isON()) { while(keyb_isON()); return; }

        if(speedup == 0) { if(keyb_isKeyPressed(0,6)) speedup = 1; }
        if(speedup) i+=0.1;
        else i+=0.01;
    }
    /*
    hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_5);
    hpg_fill_rect(0,0,24,5);
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text("Done!",0,0);
    */
    char str[40];
    sprintf(str,"Scale: 1E%d",(int)log10(axes_scale));
    hpg_set_color(hpg_stdscreen,HPG_COLOR_WHITE);
    int size = 4*strlen(str);
    hpg_fill_rect(0,0,size,5);
    hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_5);
    hpg_draw_line(0,6,size+1,6);
    hpg_draw_line(size+1,0,size+1,5);
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text(str,0,0);
}

void ask_continous_system(void)
{
    while(keyb_isAnyKeyPressed()) sys_LCDSynch();

    hpg_clear();
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text("What's the type of the system?",0,0);
    hpg_draw_text("1:Continous",0,16);
    hpg_draw_text("2:Discrete",0,24);

    while(1)
    {
        if(keyb_isKeyPressed(3,5)) //1
        {
            SYSTEM_IS_CONTINOUS = 1;
            while(keyb_isAnyKeyPressed()) sys_LCDSynch();
            hpg_clear();
            return;
        }
        else if(keyb_isKeyPressed(2,5)) //2
        {
            SYSTEM_IS_CONTINOUS = 0;
            while(keyb_isAnyKeyPressed()) sys_LCDSynch();
            hpg_clear();
            return;
        }
        sys_LCDSynch();
    }
}

int regulator_added = 0;

int add_regulator(void) //return 1 if changed
{
    while(keyb_isAnyKeyPressed()) sys_LCDSynch();

    if(regulator_added)
    {
        hpg_clear();
        hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
        hpg_draw_text("There is already a regulator.",0,0);
        hpg_draw_text("Remove it?",0,16);
        hpg_draw_text("1:YES",0,32);
        hpg_draw_text("2:NO",0,40);

        while(1)
        {
            if(keyb_isKeyPressed(3,5)) //1
            {
                restore_parsed_polys();
                regulator_added = 0;
                while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                hpg_clear();
                return 1;
            }
            else if(keyb_isKeyPressed(2,5) || keyb_isON()) //2
            {
                while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                hpg_clear();
                return 0;
            }

            sys_LCDSynch();
        }
    }

    hpg_clear();
    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
    hpg_draw_text("Regulator type?",0,0);
    if(SYSTEM_IS_CONTINOUS)
    {
        hpg_draw_text("1:PI   R(s)=k(s+a)/s",0,16);
        hpg_draw_text("2:PD   R(s)=k(s+b)",0,24);
        hpg_draw_text("3:PID  R(s)=k(s+a)(s+b)/s",0,32);
    }
    else
    {
        hpg_draw_text("1:PI   R(z)=k(2/T)(z-1)/(z+1)",0,16);
        hpg_draw_text("2:PD   R(z)=k(z-1)/(Tz)",0,24);
        hpg_draw_text("3:PID  R(z)=k(z+a)(z+b)/(z-1)z",0,32);
    }
    hpg_draw_text("Change 'k': Press arrows while",0,48);
    hpg_draw_text("   you are in the graph.",0,56);

    char string[50];

    while(1)
    {
        if(keyb_isKeyPressed(3,5)) //1 - PI
        {
            while(keyb_isAnyKeyPressed()) sys_LCDSynch();

            if(SYSTEM_IS_CONTINOUS)
            {
                double a;
                if(get_double("Input value for 'a':",&a))
                {
                    multiply_poly_regulator(num_coefs,&NUM_DEGREE,a);
                    multiply_poly_regulator(den_coefs,&DEN_DEGREE,0);

                    hpg_clear();
                    hpg_draw_text("PI regulator added!",0,0);
                    hpg_draw_text("R(s)=k(s+a)/s",0,16);
                    sprintf(string,"a=%E",a);
                    hpg_draw_text(string,0,32);

                    regulator_added = 1;

                    while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                    while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                }
            }
            else
            {
                double T;
                if(get_double("Input value for 'T':",&T))
                {
                    multiply_poly_const(num_coefs,DEN_DEGREE,2.0/T);
                    multiply_poly_regulator(num_coefs,&NUM_DEGREE,-1);
                    multiply_poly_regulator(den_coefs,&DEN_DEGREE,1);

                    hpg_clear();
                    hpg_draw_text("PI regulator added!",0,0);
                    hpg_draw_text("R(z)=k(2/T)(z-1)/(z+1)",0,16);
                    sprintf(string,"T=%E",T);
                    hpg_draw_text(string,0,32);

                    regulator_added = 1;

                    while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                    while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                }
            }
            break;
        }
        else if(keyb_isKeyPressed(2,5)) //2 - PD
        {
            while(keyb_isAnyKeyPressed()) sys_LCDSynch();

            if(SYSTEM_IS_CONTINOUS)
            {
                double b;
                if(get_double("Input value for 'b':",&b))
                {
                    multiply_poly_regulator(num_coefs,&NUM_DEGREE,b);

                    hpg_clear();
                    hpg_draw_text("PD regulator added!",0,0);
                    hpg_draw_text("R(s)=k(s+b)",0,16);
                    sprintf(string,"b=%E",b);
                    hpg_draw_text(string,0,32);

                    regulator_added = 1;

                    while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                    while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                }
            }
            else
            {
                double T;
                if(get_double("Input value for 'T':",&T))
                {
                    multiply_poly_regulator(num_coefs,&NUM_DEGREE,-1);
                    multiply_poly_const(den_coefs,DEN_DEGREE,T);
                    multiply_poly_regulator(den_coefs,&DEN_DEGREE,0);

                    hpg_clear();
                    hpg_draw_text("PD regulator added!",0,0);
                    hpg_draw_text("R(z)=k(z-1)/(Tz)",0,16);
                    sprintf(string,"T=%E",T);
                    hpg_draw_text(string,0,32);

                    regulator_added = 1;

                    while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                    while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                }
            }
            break;
        }
        else if(keyb_isKeyPressed(1,5)) //3 - PID
        {
            while(keyb_isAnyKeyPressed()) sys_LCDSynch();

            if(SYSTEM_IS_CONTINOUS)
            {
                double a, b;
                if(get_double("Input value for 'a':",&a))
                {
                    if(get_double("Input value for 'b':",&b))
                    {
                        multiply_poly_regulator(num_coefs,&NUM_DEGREE,a);
                        multiply_poly_regulator(num_coefs,&NUM_DEGREE,b);
                        multiply_poly_regulator(den_coefs,&DEN_DEGREE,0);

                        hpg_clear();
                        hpg_draw_text("PID regulator added!",0,0);
                        hpg_draw_text("R(s)=k(s+a)(s+b)/s",0,16);
                        sprintf(string,"a=%E",a);
                        hpg_draw_text(string,0,32);
                        sprintf(string,"b=%E",b);
                        hpg_draw_text(string,0,40);

                        regulator_added = 1;

                        while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                        while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                    }
                }
            }
            else
            {
                double a, b;
                if(get_double("Input value for 'a':",&a))
                {
                    if(get_double("Input value for 'b':",&b))
                    {
                        multiply_poly_regulator(num_coefs,&NUM_DEGREE,a);
                        multiply_poly_regulator(num_coefs,&NUM_DEGREE,b);
                        multiply_poly_regulator(den_coefs,&DEN_DEGREE,0);
                        multiply_poly_regulator(den_coefs,&DEN_DEGREE,-1);

                        hpg_clear();
                        hpg_draw_text("PID regulator added!",0,0);
                        hpg_draw_text("R(z)=k(z+a)(z+b)/(z-1)z",0,16);
                        sprintf(string,"a=%E",a);
                        hpg_draw_text(string,0,32);
                        sprintf(string,"b=%E",b);
                        hpg_draw_text(string,0,40);

                        regulator_added = 1;

                        while(!keyb_isAnyKeyPressed()) sys_LCDSynch();
                        while(keyb_isAnyKeyPressed()) sys_LCDSynch();
                    }
                }
            }
            break;
        }
        else if(keyb_isON()) { while(keyb_isAnyKeyPressed()) sys_LCDSynch(); break; }

        sys_LCDSynch();
    }

    hpg_clear();

    return regulator_added;
}

void draw_rlocus(void)
{
    hpg_set_indicator(HPG_INDICATOR_WAIT,HPG_COLOR_BLACK);

    sys_slowOff(); //SPEED UP!!

    auto_scale();

    draw_axes();
    draw_paths();
    draw_singularities();

    sys_slowOn(); //Normal speed

    hpg_set_indicator(HPG_INDICATOR_WAIT,HPG_COLOR_WHITE);
}

void print_current_poles(void)
{
    while(keyb_isKeyPressed(0,6)) sys_LCDSynch();

    hpg_clear();

    int i;
    double wr[MAX_DEGREE],wi[MAX_DEGREE];
    char str[100];

    int numr = solve_poly(COEFS,DEGREE,wr,wi);
    int line = 0;
    hpg_draw_text("Poles",0,0); hpg_draw_text("-----",0,6);
    for(i=0; i<numr; i++)
    {
        sprintf(str,"%f + j %f\n",wr[i],wi[i]);
        hpg_draw_text(str,0,12+6*line);
        line++;
        if(line == 11)
        {
            line = 0;
            hpg_draw_text("(cont)",107,74);
            while(!keyb_isAnyKeyPressed());
            while(keyb_isAnyKeyPressed());
            hpg_clear();
            hpg_draw_text("Poles",0,0); hpg_draw_text("-----",0,6);
        }
    }

    while(!keyb_isAnyKeyPressed());
    while(keyb_isAnyKeyPressed());
    hpg_clear();
}

void trace(void)
{
    hpg_t * rlocus = hpg_alloc_gray16_image(SCR_X,SCR_Y);
    hpg_t * rlocus_original = hpg_alloc_gray16_image(SCR_X,SCR_Y);
    hpg_t * tmp = hpg_alloc_gray16_image(SCR_X,SCR_Y);

    if(rlocus == NULL || rlocus_original == NULL || tmp == NULL)
    {
        error("ERROR: Not enough memory");
        return;
    }

    hpg_blit(hpg_stdscreen,0,0,SCR_X,SCR_Y,  rlocus,0,0);
    hpg_blit(hpg_stdscreen,0,0,SCR_X,SCR_Y,  rlocus_original,0,0);

    double k = 50.00;
    int isdrawn = 0;
    while(1)
    {
        if(keyb_isON()) break;

        if(keyb_isAnyKeyPressed())
        {
            if(keyb_isUp()||keyb_isDown()||keyb_isLeft()||keyb_isRight())
            {
                hpg_set_mode(hpg_stdscreen,HPG_MODE_XOR);
                hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);

                double val;
                double wr[MAX_DEGREE],wi[MAX_DEGREE];
                double wr_old[MAX_DEGREE],wi_old[MAX_DEGREE];
                int numr_old = 0;
                int i;

                if(isdrawn) //clear last poles
                {
                    val = get_x(k);
                    multiply_poly_add(val);
                    numr_old = solve_poly(COEFS,DEGREE,wr_old,wi_old);
                }

                if(keyb_isUp()) k+=1;
                if(keyb_isDown()) k-=1;

                if(keyb_isLeft()) k-=0.1;
                if(keyb_isRight()) k+=0.1;

                if(k < 0.0) k = 0.0;
                if(k > 99.9) k = 99.9;

                val = get_x(k);
                multiply_poly_add(val);

                int unstable = 0;

                int numr = solve_poly(COEFS,DEGREE,wr,wi);
                for(i = 0; i < numr; i++)
                {
                    plot_add(wr[i],wi[i],0,0);
                    plot_add(wr[i],wi[i],0,1);
                    plot_add(wr[i],wi[i],0,2);
                    plot_add(wr[i],wi[i],1,0);
                    plot_add(wr[i],wi[i],2,0);
                    plot_add(wr[i],wi[i],0,-1);
                    plot_add(wr[i],wi[i],0,-2);
                    plot_add(wr[i],wi[i],-1,0);
                    plot_add(wr[i],wi[i],-2,0);

                    if(SYSTEM_IS_CONTINOUS)
                    {
                        if(wr[i] >= 0.0) unstable = 1;
                    }
                    else
                    {
                        if( ((wr[i]*wr[i])+(wi[i]*wi[i])) >= 1.0) unstable = 1; //1.0^2 = 1.0
                    }
                }

                if(isdrawn) //clear last poles
                {
                    for(i = 0; i < numr_old; i++)
                    {
                        plot_add(wr_old[i],wi_old[i],0,0);
                        plot_add(wr_old[i],wi_old[i],0,1);
                        plot_add(wr_old[i],wi_old[i],0,2);
                        plot_add(wr_old[i],wi_old[i],1,0);
                        plot_add(wr_old[i],wi_old[i],2,0);
                        plot_add(wr_old[i],wi_old[i],0,-1);
                        plot_add(wr_old[i],wi_old[i],0,-2);
                        plot_add(wr_old[i],wi_old[i],-1,0);
                        plot_add(wr_old[i],wi_old[i],-2,0);
                    }
                }

                hpg_set_mode(hpg_stdscreen,HPG_MODE_PAINT);
                char str[40];
                sprintf(str,"%c K=%f",unstable?'U':'S',val);
                hpg_set_color(hpg_stdscreen,HPG_COLOR_WHITE);
                hpg_fill_rect(0,0,68,5);
                hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_5);
                hpg_draw_line(0,6,68,6);
                hpg_draw_line(6,0,6,6);
                hpg_draw_line(68,0,68,6);
                hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
                hpg_draw_text(str,0,0);

                isdrawn = 1;
            }
            else if(keyb_isKeyPressed(1,6)) //space - set k
            {
                double k_;
                if(get_double("Input value for 'k':",&k_))
                {
                    hpg_blit(rlocus,0,0,SCR_X,SCR_Y,hpg_stdscreen,0,0);

                    hpg_set_mode(hpg_stdscreen,HPG_MODE_XOR);
                    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);

                    double wr[MAX_DEGREE],wi[MAX_DEGREE];
                    int i;

                    multiply_poly_add(k_);

                    int unstable = 0;

                    int numr = solve_poly(COEFS,DEGREE,wr,wi);
                    for(i = 0; i < numr; i++)
                    {
                        plot_add(wr[i],wi[i],0,0);
                        plot_add(wr[i],wi[i],0,1);
                        plot_add(wr[i],wi[i],0,2);
                        plot_add(wr[i],wi[i],1,0);
                        plot_add(wr[i],wi[i],2,0);
                        plot_add(wr[i],wi[i],0,-1);
                        plot_add(wr[i],wi[i],0,-2);
                        plot_add(wr[i],wi[i],-1,0);
                        plot_add(wr[i],wi[i],-2,0);

                        if(SYSTEM_IS_CONTINOUS)
                        {
                            if(wr[i] >= 0.0) unstable = 1;
                        }
                        else
                        {
                            if( ((wr[i]*wr[i])+(wi[i]*wi[i])) >= 1.0) unstable = 1; //1.0^2 = 1.0
                        }
                    }

                    hpg_set_mode(hpg_stdscreen,HPG_MODE_PAINT);
                    char str[100];
                    sprintf(str,"%c K=%f",unstable?'U':'S',k_);
                    hpg_set_color(hpg_stdscreen,HPG_COLOR_WHITE);
                    hpg_fill_rect(0,0,SCR_X-1,5);
                    hpg_set_color(hpg_stdscreen,HPG_COLOR_GRAY_5);
                    hpg_draw_line(0,6,SCR_X-1,6);
                    hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
                    hpg_draw_text(str,0,0);

                    while(!keyb_isAnyKeyPressed());
                    while(keyb_isAnyKeyPressed());

                    print_current_poles();
                }

                isdrawn = 0;

                hpg_blit(rlocus,0,0,SCR_X,SCR_Y,hpg_stdscreen,0,0);
            }
            else if(keyb_isKeyPressed(0,0)) //backspace - clear
            {
                hpg_blit(rlocus,0,0,SCR_X,SCR_Y,hpg_stdscreen,0,0);

                isdrawn = 0;
            }
            else if(keyb_isAlpha()) // add regulator - clear
            {
                if(add_regulator())
                {
                    if(regulator_added)
                    {
                        draw_rlocus();
                        hpg_blit(hpg_stdscreen,0,0,SCR_X,SCR_Y,  rlocus,0,0); //save new root locus
                    }
                    else
                    {
                        auto_scale();
                        hpg_blit(rlocus_original,0,0,SCR_X,SCR_Y,hpg_stdscreen,0,0); //restore original
                        hpg_blit(rlocus_original,0,0,SCR_X,SCR_Y,  rlocus,0,0);
                    }
                }
                else
                {
                    hpg_blit(rlocus,0,0,SCR_X,SCR_Y,hpg_stdscreen,0,0);
                }
                isdrawn = 0;
            }
            else if(keyb_isKeyPressed(0,6)) //ENTER
            {
                if(isdrawn)
                {
                    hpg_blit(hpg_stdscreen,0,0,SCR_X,SCR_Y,  tmp,0,0);

                    print_current_poles();

                    hpg_blit(tmp,0,0,SCR_X,SCR_Y,  hpg_stdscreen,0,0);
                }
            }
        }

        sys_LCDSynch(); //wait for LCD Refresh
    }

    hpg_free_image(rlocus);
    hpg_free_image(rlocus_original);
    hpg_free_image(tmp);
}

void print_coefs(void)
{
    hpg_clear();

    int i;
    char str[100];
    int line = 0;
    hpg_draw_text("Numerator",0,0); hpg_draw_text("---------",0,6);
    for(i=0; i<=NUM_DEGREE; i++)
    {
        if(SYSTEM_IS_CONTINOUS)
            sprintf(str,"s^%d %f\n",NUM_DEGREE - i, num_coefs[i]);
        else
            sprintf(str,"z^%d %f\n",NUM_DEGREE - i, num_coefs[i]);
        hpg_draw_text(str,0,12+6*line);
        line++;
        if(line == 11)
        {
            line = 0;
            hpg_draw_text("(cont)",107,74);
            while(!keyb_isAnyKeyPressed());
            while(keyb_isAnyKeyPressed());
            hpg_clear();
            hpg_draw_text("Numerator",0,0); hpg_draw_text("---------",0,6);
        }
    }

    while(!keyb_isAnyKeyPressed());
    if(keyb_isON()) { hpg_clear(); return; }
    while(keyb_isAnyKeyPressed());
    hpg_clear();

    line = 0;
    hpg_draw_text("Denominator",0,0); hpg_draw_text("-----------",0,6);
    for(i=0; i<=DEN_DEGREE; i++)
    {
        if(SYSTEM_IS_CONTINOUS)
            sprintf(str,"s^%d %f\n",DEN_DEGREE - i, den_coefs[i]);
        else
            sprintf(str,"z^%d %f\n",DEN_DEGREE - i, den_coefs[i]);
        hpg_draw_text(str,0,12+6*line);
        line++;
        if(line == 11)
        {
            line = 0;
            hpg_draw_text("(cont)",107,74);
            while(!keyb_isAnyKeyPressed());
            while(keyb_isAnyKeyPressed());
            hpg_clear();
            hpg_draw_text("Denominator",0,0); hpg_draw_text("-----------",0,6);
        }
    }

    while(!keyb_isAnyKeyPressed());
    if(keyb_isON()) { hpg_clear(); return; }
    while(keyb_isAnyKeyPressed());
    hpg_clear();
}

void print_zeros_roots(void)
{
    hpg_clear();

    int i;
    double wr[MAX_DEGREE],wi[MAX_DEGREE];
    int numr;
    char str[100];
    int line;

    numr = solve_poly(num_coefs,NUM_DEGREE,wr,wi);
    line = 0;
    hpg_draw_text("Zeros",0,0); hpg_draw_text("-----",0,6);
    for(i=0; i<numr; i++)
    {
        sprintf(str,"%f + j %f\n",wr[i],wi[i]);
        hpg_draw_text(str,0,12+6*line);
        line++;
        if(line == 11)
        {
            line = 0;
            hpg_draw_text("(cont)",107,74);
            while(!keyb_isAnyKeyPressed());
            while(keyb_isAnyKeyPressed());
            hpg_clear();
            hpg_draw_text("Zeros",0,0); hpg_draw_text("-----",0,6);
        }
    }

    while(!keyb_isAnyKeyPressed());
    if(keyb_isON()) { hpg_clear(); return; }
    while(keyb_isAnyKeyPressed());
    hpg_clear();

    numr = solve_poly(den_coefs,DEN_DEGREE,wr,wi);
    line = 0;
    hpg_draw_text("Poles",0,0); hpg_draw_text("-----",0,6);
    for(i=0; i<numr; i++)
    {
        sprintf(str,"%f + j %f\n",wr[i],wi[i]);
        hpg_draw_text(str,0,12+6*line);
        line++;
        if(line == 11)
        {
            line = 0;
            hpg_draw_text("(cont)",107,74);
            while(!keyb_isAnyKeyPressed());
            while(keyb_isAnyKeyPressed());
            hpg_clear();
            hpg_draw_text("Poles",0,0); hpg_draw_text("-----",0,6);
        }
    }

    while(!keyb_isAnyKeyPressed());
    while(keyb_isAnyKeyPressed());
    hpg_clear();
}

//void LCDOff();
//void LCDOn();

int main()
{
    hpg_set_mode_gray16(0);
    hpg_clear();

    if(parse_poly() == 0)
    {
        hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);
        hpg_draw_text("Root Locus 0.3",37,5);
        hpg_draw_text("--------------",37,10);
        hpg_draw_text("by Antonio Nino Diaz",25,20);

        while(!keyb_isAnyKeyPressed());
        if(!keyb_isON())
        {
            while(keyb_isAnyKeyPressed());
        }

        hpg_clear();

        hpg_set_color(hpg_stdscreen,HPG_COLOR_BLACK);

        ask_continous_system();

        if(!keyb_isON())
        {
            print_coefs();
        }

        if(!keyb_isON())
        {
            print_zeros_roots();
        }

        while(keyb_isON());

        draw_rlocus();

        trace();
    }

    while(!keyb_isAnyKeyPressed());
    while(keyb_isAnyKeyPressed());

    hpg_cleanup();

    return(0);
}

