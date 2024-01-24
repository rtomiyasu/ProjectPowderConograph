// reflection_conditions.cc
# include <assert.h>
# include "reflection_conditions.hh"
# include "../bravais_type/BravaisType.hh"
# include "../utility_func/zstring.hh"


bool is_not_extinct_none(const Int4& h, const Int4& k, const Int4& l)
{
    return true;
}

bool is_not_extinct_4h00(const Int4& h, const Int4& k, const Int4& l)
{
	//h00: h=4n, 0k0: k=4n, 00l: l=4n
	if ( h % 4 != 0 && k == 0 && l == 0 ) return false;
	if ( h == 0 && k % 4 != 0 && l == 0 ) return false;
	if ( h == 0 && k == 0 && l % 4 != 0 ) return false;
    return true;
}

bool is_not_extinct_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//h00 : h = 2n,    0k0 : k = 2n,    00l : l = 2n
	if ( h % 2 != 0 && k == 0 && l == 0 ) return false;
	if ( h == 0 && k % 2 != 0 && l == 0 ) return false;
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
    return true;
}


bool is_not_extinct_4h00_40kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl:k+l=4n,h0l:h+l=4n,hk0:h+k=4n,h00:h=4n,0k0:k=4n,00l:l=4n
	if ( !is_not_extinct_4h00(h, k, l) ) return false;
	if ( h == 0 && (k + l) % 4 != 0 ) return false;
	if ( k == 0 && (h + l) % 4 != 0 ) return false;
	if ( l == 0 && (h + k) % 4 != 0 ) return false;
    return true;
}

bool is_not_extinct_2h00_20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 2n,    h0l : h+l = 2n,    hk0 : h+k = 2n,    h00 : h = 2n,      0k0 : k = 2n,      00l : l = 2n
	if ( !is_not_extinct_2h00(h, k, l) ) return false;
	if ( h == 0 && (k + l) % 2 != 0 ) return false;
	if ( k == 0 && (h + l) % 2 != 0 ) return false;
	if ( l == 0 && (h + k) % 2 != 0 ) return false;
    return true;
}

bool is_not_extinct_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl:h,l=2n,  hkh:h,k=2n,  hkk:h,k=2n
	if ( k == h && !(h % 2 == 0 && l % 2 == 0 )) return false;
	if ( l == h && !(h % 2 == 0 && k % 2 == 0 ) ) return false;
	if ( l == k && !(h % 2 == 0 && k % 2 == 0 )    ) return false;
    return true;
}

bool is_not_extinct_20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl:k,l=2n,  h0l:h,l=2n,  hk0:h,k=2n
	if ( h == 0 && !(k % 2 == 0 && l % 2 == 0 )) return false;
	if ( k == 0 && !(h % 2 == 0 && l % 2 == 0 )) return false;
	if ( l == 0 && !(h % 2 == 0 && k % 2 == 0 )) return false;
    return true;
}
bool is_not_extinct_40kl_2hhl_4h00(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl:k+l=4n,  h0l:h+l=4n,  hk0:h+k=4n,  hhl:h,l=2n,  hkh:h,k=2n,  hkk:h,k=2n,  h00:h=4n,  0k0:k=4n,  00l:l=4n
	if ( !is_not_extinct_4h00_40kl(h, k, l) ) return false;
	if ( !is_not_extinct_2hhl(h, k, l) ) return false;
    return true;
}
bool is_not_extinct_4hhl_4h00(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : 2h+l = 4n,    hkh : 2h+k = 4n,    hkk : h+2k = 4n,    h00 : h = 4n,      0k0 : k = 4n,      00l : l = 4n
	if ( !is_not_extinct_4h00(h, k, l) ) return false;
	if ( k == h && (2*h+l) % 4 !=0 ) return false;
	if ( l == h && (2*h+k) % 4 !=0 ) return false;
	if ( l == k && (2*k+h) % 4 !=0 ) return false;
    return true;
}

bool is_not_extinct_2hhl_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : l = 2n,    hkh : k = 2n,    hkk : h = 2n,    h00 : h = 2n,      0k0 : k = 2n,      00l : l = 2n
	if ( !is_not_extinct_2h00(h, k, l) ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
	if ( l == h && (k) % 2 !=0 ) return false;
	if ( l == k && (h) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_2hhl_2h00_20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl:k+l=2n,   h0l:h+l=2n,   hk0:h+k=2n,    hhl:l=2n,  hkh:k=2n,  hkk:h=2n,  h00:h=2n,   0k0:k=2n,  00l:l=2n
	if ( !is_not_extinct_2h00_20kl(h, k, l) ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
	if ( l == h && (k) % 2 !=0 ) return false;
	if ( l == k && (h) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_2hk0_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h = 2n      0kl : k = 2n      h0l : l = 2n     h00 : h = 2n    0k0 : k = 2n    00l : l = 2n
	if ( !is_not_extinct_2h00(h, k, l) ) return false;
	if ( l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && (k) % 2 !=0 ) return false;
	if ( k == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_2hk0mirror_2h00(const Int4& h, const Int4& k, const Int4& l)
{
        //h00 : h = 2n    0k0 : k = 2n    00l : l = 2n
        if ( !is_not_extinct_2h00(h, k, l) ) return false;
        //hk0 : k = 2n      0kl : l = 2n      h0l : h = 2n
        if ( l == 0 && (k) % 2 !=0 ) return false;
        if ( h == 0 && (l) % 2 !=0 ) return false;
        if ( k == 0 && (h) % 2 !=0 ) return false;
    return true;
}
//Trigonal: tr
bool is_not_extinct_tr_3000l(const Int4& h, const Int4& k, const Int4& l)
{
//000l : l = 3n
	if ( h == 0 && k == 0 && (l) % 3 !=0 ) return false;
    return true;
}
bool is_not_extinct_tr_2hmh0l(const Int4& h, const Int4& k, const Int4& l)
{
//h-h0l : l = 2n    h0-hl : l = 2n    0h-hl : l = 2n
	if ( k == -h && (l) % 2 !=0 ) return false;
	if ( k == 0 && (l) % 2 !=0 ) return false;
	if ( h == 0 && (l) % 2 !=0 ) return false;
    return true;
}
//Hexagonal : hex
bool is_not_extinct_hex_6000l(const Int4& h, const Int4& k, const Int4& l)
{
//000l : l = 6n
	if ( h == 0 && k == 0 && (l) % 6 !=0 ) return false;
    return true;
}
bool is_not_extinct_hex_3000l(const Int4& h, const Int4& k, const Int4& l)
{
//000l : l = 3n
	if ( h == 0 && k == 0 && (l) % 3 !=0 ) return false;
    return true;
}
bool is_not_extinct_hex_2000l(const Int4& h, const Int4& k, const Int4& l)
{
//000l : l = 2n
	if ( h == 0 && k == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_hex_2hhm2hl_2hmh0l_2000l(const Int4& h, const Int4& k, const Int4& l)
{
//hh-2hl : l = 2n,    h-2hhl : l = 2n,    -2hhhl : l = 2n,    h-h0l : l = 2n,    h0-hl : l = 2n,    0h-hl : l = 2n
	if ( k == h && (l) % 2 !=0 ) return false;
	if ( k == (-2 * h) && (l) % 2 !=0 ) return false;
	if ( h == (-2 * k) && (l) % 2 !=0 ) return false;
	if ( k == (-h) && (l) % 2 !=0 ) return false;
	if ( k == 0 && (l) % 2 !=0 ) return false;
	if ( h == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_hex_2hmh0l(const Int4& h, const Int4& k, const Int4& l)
{
//h-h0l : l = 2n,    h0-hl : l = 2n,    0h-hl : l = 2n
	if ( k == (-h) && (l) % 2 !=0 ) return false;
	if ( k == 0 && (l) % 2 !=0 ) return false;
	if ( h == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_2hhm2hl(const Int4& h, const Int4& k, const Int4& l)
{
//hh-2hl : l = 2n,    h-2hhl : l = 2n,    -2hhhl : l = 2n
	if ( k == h && (l) % 2 !=0 ) return false;
	if ( k == (-2 * h) && (l) % 2 !=0 ) return false;
	if ( h == (-2 * k) && (l) % 2 !=0 ) return false;
    return true;
}
//Tetragonal
bool is_not_extinct_400l(const Int4& h, const Int4& k, const Int4& l)
{
	//00l:l=4n
	if ( h == 0 && k == 0 && (l) % 4 !=0 ) return false;
    return true;
}

bool is_not_extinct_200l(const Int4& h, const Int4& k, const Int4& l)
{
	//00l:l=2n
	if ( h == 0 && k == 0 && (l) % 2 !=0 ) return false;
    return true;
}

bool is_not_extinct_2hk0_200l(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,  00l : l = 2n
	if ( !is_not_extinct_200l(h, k, l) ) return false;
	if ( l == 0 && (h + k) %  2 != 0) return false;
    return true;
}

bool is_not_extinct_2h00_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	//h00 : h = 2n,    0k0 : k = 2n
	if ( h % 2 != 0 && k == 0 && l == 0 ) return false;
	if ( h == 0 && k % 2 != 0 && l == 0 ) return false;
    return true;
}

bool is_not_extinct_400l_2h00_0k0(const Int4& h, const Int4& k, const Int4& l)
{
	//00l : l = 4n,    h00 : h = 2n,    0k0 : k = 2n
	if ( !is_not_extinct_2h00_20k0(h, k, l) ) return false;
	if ( h == 0 && k == 0 && (l) % 4 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : l = 2n,        0kl : l = 2n
	if ( k == 0 && (l) % 2 !=0 ) return false;
	if ( h == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2h0l_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : l = 2n,        0kl : l = 2n,     h00 : h = 2n,    0k0 : k = 2n
	if ( !is_not_extinct_t_2h0l(h, k, l) ) return false;
	if ( k == 0 && l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && l == 0 && (k) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2ph0l_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h+l = 2n, 0kl : k+l = 2n ,h00 : h = 2n,    0k0 : k = 2n
	if ( k == 0 && (h + l) % 2 !=0 ) return false;
	if ( h == 0 && (k + l) % 2 !=0 ) return false;
	if ( k == 0 && l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && l == 0 && (k) % 2 !=0 ) return false;
    return true;
}

bool is_not_extinct_t_2h0l_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : l = 2n,     0kl : l = 2n,     hhl : l = 2
	if ( !is_not_extinct_t_2h0l(h, k, l) ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : l = 2n
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}

bool is_not_extinct_t_2h00_2h0l_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h = 2n,    0kl : k = 2n,    hhl : l = 2n,     h00 : h = 2n,     0k0 : k = 2n
	if ( k == 0 && l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && l == 0 && (k) % 2 !=0 ) return false;
	if ( k == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && (k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}

bool is_not_extinct_t_2h00_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : l = 2n,     h00 : h = 2n,     0k0 : k = 2n
	if ( k == 0 && l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && l == 0 && (k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2phk0_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,  h0l : h+l = 2n,  0kl : k+l = 2n,  hhl : l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == 0 && (h + l) % 2 !=0 ) return false;
	if ( h == 0 && (l + k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2phk0_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,   h0l : h = 2n,    0kl : k = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && (k) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2ph0l_2hhl_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h+l = 2n,   0kl : k+l = 2n,  hhl : l = 2n,     h00 : h = 2n,     0k0 : k = 2n
	if ( k == 0 && (h + l) % 2 !=0 ) return false;
	if ( h == 0 && (l + k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
	if ( k == 0 && l == 0 && (h) % 2 !=0 ) return false;
	if ( h == 0 && l == 0 && (k) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2phk0(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2phk0_2h0l_2hll(const Int4& h, const Int4& k, const Int4& l)
{
    //hk0 : h+k = 2n,    h0l : l = 2n,     0kl : l = 2n,     hhl : l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == 0 && (l ) % 2 !=0 ) return false;
	if ( h == 0 && (l ) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2phk0_2hll(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,    hhl : l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_23phk0(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,  h0l : h+l = 2n,  0kl : k+l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == 0 && (h + l) % 2 !=0 ) return false;
	if ( h == 0 && (l + k) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_21phk0_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,  hhl : l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == h && (l) % 2 !=0 ) return false;
    return true;
}

bool is_not_extinct_t_21phk0_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h+k = 2n,    h0l : l = 2n,       0kl : l = 2n
	if ( l == 0 && (h + k) % 2 !=0 ) return false;
	if ( k == 0 && (l) % 2 !=0 ) return false;
	if ( h == 0 && (l) % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_21hk0_200l(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h,k = 2n,  00l : l = 4n
	if ( l == 0 && !(h % 2 == 0 && (k) % 2 == 0 ) ) return false;
	if ( k == 0 && h == 0 && (l) % 4 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_2ah0l(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h, l = 2n,    0kl : k, l = 2n
	if ( k == 0 && !(h % 2 == 0 && (l) % 2 == 0 ) ) return false;
	if ( h == 0 && !(k % 2 == 0 && (l) % 2 == 0 ) ) return false;
    return true;
}
bool is_not_extinct_t_4phhl_2hmh0(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : 2h+l = 4n,  h-h0 : h = 2n
	if ( k == h && (2 * h + l) % 4 !=0 ) return false;
	if ( l == 0 && k == -h && h % 2 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_4phhl_2hmh0_2ah0l(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : 2h+l = 4n,  h-h0 : h = 2n,    h0l : h, l = 2n,    0kl : k, l = 2n
	if ( k == h && (2 * h + l) % 4 !=0 ) return false;
	if ( l == 0 && k == -h && h % 2 !=0 ) return false;
	if ( k == 0 && !(h % 2 == 0 && (l) % 2 == 0 ) ) return false;
	if ( h == 0 && !(k % 2 == 0 && (l) % 2 == 0 ) ) return false;
    return true;
}
bool is_not_extinct_t_2hk0_4phhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h,k = 2n,  hhl : 2h+l = 4n
	if ( l == 0 && !(h % 2 == 0 && (k) % 2 == 0 ) ) return false;
	if ( k == h && (2 * h + l) % 4 !=0 ) return false;
    return true;
}
bool is_not_extinct_t_23ahk0_4phhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h,k = 2n,  h0l : h, l = 2n,    0kl : k, l = 2n,  hhl : 2h+l = 4n
	if ( l == 0 && !(h % 2 == 0 && (k) % 2 == 0 ) ) return false;
	if ( k == 0 && !(h % 2 == 0 && (l) % 2 == 0 ) ) return false;
	if ( h == 0 && !(l % 2 == 0 && (k) % 2 == 0 ) ) return false;
	if ( k == h && (2 * h + l) % 4 !=0 ) return false;
    return true;
}
bool is_not_extinct_rhom_2hhl(const Int4& h, const Int4& k, const Int4& l)
{
	//hhl : l = 2n,    hkh: k = 2n,   hkk : h = 2n
	if ( k == h && l % 2 != 0 ) return false;
	if ( l == h && k % 2 != 0 ) return false;
	if ( l == k && h % 2 != 0 ) return false;
    return true;
}

bool is_not_extinct_rhom_hex_2hmh0l(const Int4& h, const Int4& k, const Int4& l)
{
	// h-h0l : l = 2n,    h0-hl : l = 2n,    0h-hl : l = 2n
	if ( k == - h && l % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( h == 0 && l % 2 != 0 ) return false;
    return true;
}

//Monoclinic P
bool is_not_extinct_mono_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	// h00:h=2n
	if ( k == 0 && l == 0 && h % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	// 0k0:k=2n
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_200l(const Int4& h, const Int4& k, const Int4& l)
{
	// 00l:l=2n
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_20kl(const Int4& h, const Int4& k, const Int4& l)
{
	// 0kl:k=2n
	if ( h == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	// h0l:l=2n
	if ( k == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2hk0(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:h=2n
	if ( l == 0 && h % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2l0kl(const Int4& h, const Int4& k, const Int4& l)
{
	// 0kl:l=2n
	if ( h == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2hh0l(const Int4& h, const Int4& k, const Int4& l)
{
	// h0l:h=2n
	if ( k == 0 && h % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2khk0(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:k=2n
	if ( l == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2hk0_200l(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:h=2n,   00l : l = 2n
	if ( l == 0 && h % 2 != 0 ) return false;
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2phk0_200l(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:h+k=2n,   00l : l = 2n
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2dhk0_200l(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:k=2n,   00l : l = 2n
	if ( l == 0 && k % 2 != 0 ) return false;
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2p0kl(const Int4& h, const Int4& k, const Int4& l)
{
	// 0kl:k+l=2n
	if ( h == 0 && (k + l) % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2ph0l(const Int4& h, const Int4& k, const Int4& l)
{
	// h0l:h+l=2n
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2phk0(const Int4& h, const Int4& k, const Int4& l)
{
	// hk0:h+k=2n
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2h0l_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : l = 2n,      0k0 : k = 2n
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2ph0l_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h + l = 2n,      0k0 : k = 2n
	if ( k == 0 && ( h +l ) % 2 != 0 ) return false;
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2dh0l_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h = 2n,      0k0 : k = 2n
	if ( k == 0 && h % 2 != 0 ) return false;
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_20kl_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k = 2n,      h00 : h = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l == 0 && h % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2d0kl_2h00(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2n,      h00 : h = 2n
	if ( h == 0 && l % 2 != 0 ) return false;
	if ( k == 0 && l == 0 && h % 2 != 0 ) return false;
    return true;
}
bool is_not_extinct_mono_2p0kl_20k0(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k + l = 2n,      0k0 : k = 2n
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
    return true;
}
//Orthorhombic
bool is_not_extinct_orth_2dh00(const Int4& h, const Int4& k, const Int4& l)
{
//h00 : h = 2n,    0k0 : k = 2n,     00l : l = 2n
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l == 0 && h % 2 != 0 ) return false;
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool is_not_extinct_orth_23phk0(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0:h+k=2n,h0l:h+l=2n,0kl:k+l=2n
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	return true;
}
bool is_not_extinct_orth_2d30kl(const Int4& h, const Int4& k, const Int4& l)
{
//0kl : k = 2n,縲�    h0l : l = 2n,  縲�  hk0 : h = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( l == 0 && h  % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl:l = 2n,  h0l:l=2n, hk0:h+k=2n
	if ( h == 0 && l % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l:l=2n
	if ( k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_200l(const Int4& h, const Int4& k, const Int4& l)
{
	//00l:l = 2n
	if ( h == 0 && k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_22h00(const Int4& h, const Int4& k, const Int4& l)
{
	// h00 : k = 2n,    0k0 : k = 2n.
	if ( h == 0 && l == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_220kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2n,       h0l : l = 2n
	if ( h == 0 && l % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_2d20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k = 2n,       h0l : h = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_2p20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 2n    h0l : h+l = 2n
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_220kl_2phk0(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k = 2n,   縲�縲� h0l : h = 2n,  縲�縲�  hk0 : h+k = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && h % 2 != 0 ) return false;
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_22d0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2n,       h0l : h = 2n
	if ( h == 0 && l % 2 != 0 ) return false;
	if ( k == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_22dd0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k = 2n,       h0l : l = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_2p0kl_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl: k+l=2n ,   h0l: l=2n
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_2ph0l_2hk0(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l: h+l=2n ,   hk0: h=2n
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
	if ( l == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_2ph0l(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l: h+l=2n
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_2p0kl_2dh0l(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 2n,    h0l : h = 2n
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	if ( k == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_22p0kl_2hk0(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 2n ,縲�縲�  h0l : h+l = 2n,    hk0 : h = 2n
	if ( h == 0 && ( k + l ) % 2 != 0 ) return false;
	if ( k == 0 && ( h + l ) % 2 != 0 ) return false;
	if ( l == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_230kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2,    縲�縲� h0l : l = 2n ,   縲�hk0 : h = 2n
	if ( h == 0 && l % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( l == 0 && h % 2 != 0 ) return false;
	return true;
}
bool standard_function2_for_abc_220kl_2phk0(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k = 2n,    縲�h0l : l = 2n,   縲� hk0 : h+k = 2n
	if ( h == 0 && k % 2 != 0 ) return false;
	if ( k == 0 && l % 2 != 0 ) return false;
	if ( l == 0 && ( h + k ) % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_210kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2n
	if ( h == 0 && l % 2 != 0 ) return false;
	return true;
}
bool standard_function_for_abc_2ahk0(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h,k = 2n
	if ( l == 0 && !( h % 2 == 0 && k % 2 == 0 ) ) return false;
	return true;
}
bool standard_function_for_abc_2ahk0_20kl(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h, k = 2n,   0kl : l = 2n
	if ( l == 0 && !( h % 2 == 0 && k % 2 == 0 ) ) return false;
	if ( h == 0 &&  l % 2 != 0  ) return false;
	return true;
}
bool standard_function_for_abc_2ahk0_2h0l(const Int4& h, const Int4& k, const Int4& l)
{
	//hk0 : h, k = 2n,   h0l : l = 2n
	if ( l == 0 && !( h % 2 == 0 && k % 2 == 0 ) ) return false;
	if ( k == 0 &&  l % 2 != 0  ) return false;
	return true;
}
bool standard_function_for_abc_220kl_2ahk0(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : l = 2n,   縲� h0l : l = 2n,    hk0 : h,k = 2n
	if ( h == 0 &&  l % 2 != 0  ) return false;
	if ( k == 0 &&  l % 2 != 0  ) return false;
	if ( l == 0 && !( h % 2 == 0 && k % 2 == 0 ) ) return false;
	return true;
}
bool standard_function_for_abc_42p0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 4n,           h0l : h+l = 4n
	if ( h == 0 &&  ( k+l) % 4 != 0  ) return false;
	if ( k == 0 &&  ( h+l) % 4 != 0  ) return false;
	return true;
}
bool standard_function_for_abc_43p0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k+l = 4n,           h0l : h+l = 4n,  hk0 : h+k = 4n
	if ( h == 0 &&  ( k+l) % 4 != 0  ) return false;
	if ( k == 0 &&  ( h+l) % 4 != 0  ) return false;
	if ( l == 0 &&  ( h+k) % 4 != 0  ) return false;
	return true;
}
bool standard_function_for_abc_2a0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k,l = 2n     h0l : h,l = 2n
	if ( h == 0 &&  !( k % 2 == 0 && l % 2 == 0 )  ) return false;
	if ( k == 0 &&  !( h % 2 == 0 && l % 2 == 0 )  ) return false;
	return true;
}
bool standard_function_for_abc_2ah0l(const Int4& h, const Int4& k, const Int4& l)
{
	//h0l : h,l = 2n
	if ( k == 0 &&  !( h % 2 == 0 && l % 2 == 0 )  ) return false;
	return true;
}
bool standard_function_for_abc_23a0kl(const Int4& h, const Int4& k, const Int4& l)
{
	//0kl : k, l = 2n  h0l : h, l = 2n  hk0 : h,k = 2n
	if ( h == 0 &&  !( k % 2 == 0 && l % 2 == 0 )  ) return false;
	if ( k == 0 &&  !( h % 2 == 0 && l % 2 == 0 )  ) return false;
	if ( l == 0 &&  !( h % 2 == 0 && k % 2 == 0 )  ) return false;
	return true;
}
bool special_reflection_conditions_3h_4phkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (A,face) hkl : h=2n+1 or h+k+l=4n
	if( (h+k+l-2) % 4 == 0 ) return false;
	return true;
}
bool special_reflection_conditions_3l_4p2hl(const Int4& h, const Int4& k, const Int4& l)
{
	// (A,body) hkl : l=2n+1 or 2h+l=4n
	if( (2*h+l-2) % 4 == 0 ) return false;
	return true;
}
bool special_reflection_conditions_3h_6ahkl_4ahkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (B,face) hkl:h=2n+1 or h,k,l=4n+2 or h,k,l=4n
	if( h % 2 == 0 && !( (h - k) % 4 == 0 && (k - l) % 4 == 0 ) ) return false;
	return true;
}
bool special_reflection_conditions_3l_2ahk_4phkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (B,body) hkl : l=2n+1or  h,k=2n, h+k+l=4n
	if( l % 2 == 0 )
	{
	    if( h % 2 != 0 ) return false;
	    if( (h + k + l) % 4 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_3l_2h(const Int4& h, const Int4& k, const Int4& l)
{
	// (C) hkl : l=2n+1    or  h=2n
	if( h % 2 != 0 && l % 2 == 0 ) return false;
	return true;
}
bool special_reflection_conditions_2l_3shk_3s2hk(const Int4& h, const Int4& k, const Int4& l)
{
	// (D) hkil : l=2n  or h-k=3n+1 or h-k=3n+2
	if( (h - k) % 3 == 0 && l % 2 != 0 ) return false;
	return true;
}
bool special_reflection_conditions_oddh_oddk_3lhkil(const Int4& h, const Int4& k, const Int4& l)
{
	// (E) hkil : h=2n+1  or k=2n+1 or l=3n
	if( h % 2 == 0 && k % 2 == 0 && l % 3 != 0 ) return false;
	return true;
}
bool special_reflection_conditions_2h_2k_2hhl_2hkh_2hkk(const Int4& h, const Int4& k, const Int4& l)
{
	// (F) hkl: h=2n or k=2n or l=2n,  hhl: l=2n,  hkh: k=2n,  hkk: h=2n
	if( k == l && h % 2 != 0 ) return false;
	if( h == l && k % 2 != 0 ) return false;
	if( h == k && l % 2 != 0 ) return false;
	if( h % 2 != 0 && k % 2 != 0 && l % 2 != 0 ) return false;
	return true;
}
bool special_reflection_conditions_2h_2k_2lhkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (F) hkl: h=2n or k=2n or l=2n
	if( h % 2 != 0 && k % 2 != 0 && l % 2 != 0 ) return false;
	return true;
}
bool special_reflection_conditions_oddh_oddk_oddl_4phkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (F) hkl: h=2n+1 or k=2n+1 or l=2n+1 or h+k+l=4n
	if( h % 2 == 0 && k % 2 == 0 && l % 2 == 0 && (h+k+l) % 4 != 0) return false;
	return true;
}
bool special_reflection_conditions_2phkl_oddhkl_4hkl_6hkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (G) hkl:h+k+l=2n or either one of h, k, l is 2n+1, 4n and 4n+2
	if( (h+k+l) % 2 != 0 )
	{
	   if( h % 4 != 0 && k % 4 != 0 && l % 4 != 0 ) return false;
	   if( (h-2) % 4 != 0 && (k-2) % 4 != 0 && (l-2) % 4 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_22hkl_2oddhkl_4k_6l_6ahkl_4hkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (H) hkl: two of h,k,l are odd or either one of h,k,l is 2n+1,k=4n and l=4n+2 or h,k,l=4n+2 or h,k,l=4n
	const Int4 num_odd = (abs(h) % 2) + (abs(k) % 2) + (abs(l) % 2);
	if( num_odd < 1 )
	{
	    if( (h - k) % 4 != 0 || (k - l) % 4 != 0 ) return false;
	}
	else if( num_odd < 2 )
	{
	   if( h % 4 != 0 && k % 4 != 0 && l % 4 != 0 ) return false;
	   if( (h - 2) % 4 != 0 && (k - 2) % 4 != 0 && (l - 2) % 4 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_oddh_6ahkl_4ahkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (I) hkl:h=2n+1 or h,k,l=4n+2 or h,k,l=4n
	if( h % 2 == 0 && k % 2 == 0 && l % 2 == 0 )
	{
	    if( (h - k) % 4 != 0 || (k - l) % 4 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_2ahkl_4phkl_8h_12k_6phkl_oddahk_9ahk_4l_11ahk_4l_9h_7k_4l_11h(const Int4& h, const Int4& k, const Int4& l)
{
	// (J) hkl:h,k,l=2n,h+k+l=4n or h=8n, k=8n+4 and h+k+l=4n+2 or h,k=2n+1, l=4n+2 or h,k=8n+1, l=4n or h,k=8n+3,
	//l=4n or h=8n+1, k=8n-1, l=4n or h=8n+3, k=8n-3, l=4n
	if( h % 2 == 0 )
	{
	    if( k % 2 == 0 )  // All even
	    {
	        if( (h + k + l) % 4 != 0 )
	        { // h+k+l = 4n+2
	           if( h % 8 != 0 && k % 8 != 0 && l % 8 != 0 ) return false;
	           if( (h - 4) % 8 != 0 && (k - 4) % 8 != 0 && (l - 4) % 8 != 0 ) return false;
	        }
	    }
	    else // k, l are odd
	    {
	        if( h % 4 == 0 && (k + l) % 8 != 0 && !( (k - l) % 8 == 0 && abs(k) % 8 < 4 ) ) return false;
	    }
	}
	else if( k % 2 == 0 ) // h, l are odd
	{
	    if( k % 4 == 0 && (h + l) % 8 != 0 && !( (h - l) % 8 == 0 && abs(h) % 8 < 4 ) ) return false;
	}
	else // h, k are odd
	{
	    if( l % 4 == 0 && (h + k) % 8 != 0 && !( (h - k) % 8 == 0 && abs(h) % 8 < 4 ) ) return false;
	}

	return true;
}
bool special_reflection_conditions_2hk_4phkl_8h_12k_6phkl_oddahk_6l_9h_11k_4l_15h_13k_9h_14k_4l_15h_11k_4l(const Int4& h, const Int4& k, const Int4& l)
{
	// (J2) hkl:h,k=2n,h+k+l=4n or h=8n, k=8n+4 and h+k+l=4n+2 or h,k=2n+1, l=4n+2 or h=8n+1, k=8n+3, l=4n
	//or h=8n+7, k=8n+5, l=4n or h=8n+1, k=8n+5, l=4n or h=8n+7, k=8n+3, l=4n
	if( h % 2 == 0 )
	{
	    if( k % 2 == 0 )  // All even
	    {
	        if( (h + k + l) % 4 != 0 )
	        { // h+k+l = 4n+2
	           if( h % 8 != 0 && k % 8 != 0 && l % 8 != 0 ) return false;
	           if( (h - 4) % 8 != 0 && (k - 4) % 8 != 0 && (l - 4) % 8 != 0 ) return false;
	        }
	    }
	    else // k, l are odd
	    {
	        if( h % 4 == 0 && (k - l - 4) % 8 != 0 && (k + l - 4) % 8 == 0 ) return false;
	    }
	}
	else if( k % 2 == 0 ) // h, l are odd
	{
	    if( k % 4 == 0 && (h - l - 4) % 8 != 0 && (h + l - 4) % 8 != 0 ) return false;
	}
	else // h, k are odd
	{
	    if( l % 4 == 0 && (h - k - 4) % 8 != 0 && (h + k - 4) % 8 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_oddh_4h_oddh_4phkl_hhl(const Int4& h, const Int4& k, const Int4& l)
{
	// (K) hkl: h=2n+1  or h=4n, hhl: h=2n+1  or h+k+l=4n
	if( h % 2 == 0 && k % 2 == 0 && l % 2 == 0 )
	{
	    if( ( h == k || k == l || h == l ) && (h + k + l) % 4 != 0 ) return false;
	    if( (h - 2) % 4 == 0 && (k - 2) % 4 == 0 && (l - 2) % 4 == 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_oddh_oddk_oddl_4hkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (K) hkl: h=2n+1 or k=2n+1 or l=2n+1 or h=4n or k=4n or l=4n
	if( (h - 2) % 4 == 0 && (k - 2) % 4 == 0 && (l - 2) % 4 == 0 ) return false;
	return true;
}
bool special_reflection_conditions_4phl_4pkh_4plk_4plh_4phk_4pkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (L) hkl: 2h+l=4n or 2k+h=4n or 2l+k=4n or 2l+h=4n or 2h+k=4n or 2k+l=4n
	//(Equivalently, {u_1,u_2,u_3\} mod 4 != {0,1,1}, {0,1,3}, {0,3,3}, {2,2,2})
	if( h % 2 != 0 )
	{ // two odds and one even
	    if( k % 4 != 0 && l % 4 != 0 ) return false;
	}
	else if( k % 2 != 0 )
	{ // only h is even
	    if( h % 4 != 0 ) return false;
	}
	else
	{ // all even
	    if( (h - 2) % 4 != 0 || (k - 2) % 4 != 0 || (l - 2) % 4 != 0 ) return false;
	}
	return true;
}
bool special_reflection_conditions_2ahk_4phkl_oddahk_6l_8h_12k_6phkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (M) hkl: h,k=2n, h+k+l=4n  or h,k=2n+1, l=4n+2 or h=8n, k=8n+4 and  h+k+l=4n+2
	if( h % 2 != 0 )
	{ // two odds and one even
	    if( (k - 2) % 4 != 0 && (l - 2) % 4 != 0 ) return false;
	}
	else if( k % 2 != 0 )
	{ // only h is even
	    if( (h - 2) % 4 != 2 ) return false;
	}
	else
	{ // all even
	    if( (h + k + l) % 2 != 0 ) return false;
//	    if( (h + k + l) % 4 != 0 )
//	    {
	        if( h % 8 != 0 && k % 8 != 0 && l % 8 != 0 ) return false;
	        if( (h - 4) % 8 != 0 && (k - 4) % 8 != 0 && (l - 4) % 8 != 0 ) return false;
//	    }
	}

	return true;
}
bool special_reflection_conditions_oddahk_6l_4ahkl(const Int4& h, const Int4& k, const Int4& l)
{
	// (N) hkl: h,k=2n+1, l=4n+2  or h,k,l=4n
	if( h % 2 != 0 )
	{ // two odds and one even
	    if( (k - 2) % 4 != 0 && (l - 2) % 4 != 0 ) return false;
	}
	else if( k % 2 != 0 )
	{ // only h is even
	    if( (h - 2) % 4 != 0 ) return false;
	}
	else
	{ // all even
	    if( h % 4 != 0 || k % 4 != 0 || l % 4 != 0 ) return false;
	}
	return true;
}

static const Int4 DATA_NUM_CUBIC_F = 7;
static const Int4 DATA_NUM_CUBIC_I = 13;
static const Int4 DATA_NUM_CUBIC_P = 12;
static const Int4 DATA_NUM_HEXAGONAL = 9;
static const Int4 DATA_NUM_TETRAGONAL_P = 23;
static const Int4 DATA_NUM_TETRAGONAL_I = 11;
static const Int4 DATA_NUM_ORTHORHOMBIC_P =71 ; //95;
static const Int4 DATA_NUM_ORTHORHOMBIC_C =9;
static const Int4 DATA_NUM_ORTHORHOMBIC_F =7;
static const Int4 DATA_NUM_ORTHORHOMBIC_I =8;
static const Int4 DATA_NUM_RHOMBOHEDRAL_RHOM_AXIS = 2;
static const Int4 DATA_NUM_RHOMBOHEDRAL_HEX_AXIS = 2;
static const Int4 DATA_NUM_MONOCLINIC_P_A_AXIS =8;
static const Int4 DATA_NUM_MONOCLINIC_P_B_AXIS =8;
static const Int4 DATA_NUM_MONOCLINIC_P_C_AXIS =8;
static const Int4 DATA_NUM_MONOCLINIC_B_A_AXIS =2;
static const Int4 DATA_NUM_MONOCLINIC_B_B_AXIS =2;
static const Int4 DATA_NUM_MONOCLINIC_B_C_AXIS =2;
static const Int4 DATA_NUM_TRICLINIC =1;


Int4 putNumberOfTypesOfSystematicAbsences(const BravaisType& type)
{
    if( type.enumBravaisType() == Cubic_F )
    {
    	return DATA_NUM_CUBIC_F;
    }
    if( type.enumBravaisType() == Cubic_I )
    {
    	return DATA_NUM_CUBIC_I;
    }
    if ( type.enumBravaisType() == Cubic_P )
    {
    	return DATA_NUM_CUBIC_P;
    }
    if ( type.enumBravaisType() == Hexagonal )
    {
        return DATA_NUM_HEXAGONAL;
    }
    if ( type.enumBravaisType() == Tetragonal_P )
    {
        return DATA_NUM_HEXAGONAL;
    }
    if ( type.enumBravaisType() == Tetragonal_I )
    {
        return DATA_NUM_HEXAGONAL;
    }
    if ( type.enumBravaisType() == Rhombohedral )
    {
           if( type.enumRHaxis() == Rho_Axis )
            {
                    return DATA_NUM_RHOMBOHEDRAL_RHOM_AXIS;
            }
            if( type.enumRHaxis() == Hex_Axis )
            {
                    return DATA_NUM_RHOMBOHEDRAL_HEX_AXIS;
            }
     }
    if ( type.enumBravaisType() == Orthorhombic_P )
    {
        return DATA_NUM_ORTHORHOMBIC_P;
    }
    if ( type.enumBravaisType() == Orthorhombic_C )
    {
        return DATA_NUM_ORTHORHOMBIC_C;
    }
    if ( type.enumBravaisType() == Orthorhombic_F )
    {
        return DATA_NUM_ORTHORHOMBIC_F;
    }
    if ( type.enumBravaisType() == Orthorhombic_I )
    {
        return DATA_NUM_ORTHORHOMBIC_I;
    }    
    if ( type.enumBravaisType() == Monoclinic_P )
    {
        if( type.enumABCaxis() == A_Axis )
        {
                return DATA_NUM_MONOCLINIC_P_A_AXIS;
        }
        if( type.enumABCaxis() == B_Axis )
        {
                return DATA_NUM_MONOCLINIC_P_B_AXIS;
        }
        if( type.enumABCaxis() == C_Axis )
        {
                return DATA_NUM_MONOCLINIC_P_C_AXIS;
        }
    }
    if ( type.enumBravaisType() == Monoclinic_B )
    {
        if( type.enumABCaxis() == A_Axis )
        {
                return DATA_NUM_MONOCLINIC_B_A_AXIS;
        }
        if( type.enumABCaxis() == B_Axis )
        {
                return DATA_NUM_MONOCLINIC_B_B_AXIS;
        }
        if( type.enumABCaxis() == C_Axis )
        {
                return DATA_NUM_MONOCLINIC_B_C_AXIS;
        }
    }
    if ( type.enumBravaisType() == Triclinic )
    {
    	return DATA_NUM_TRICLINIC;
    }

        assert( false );
    return -1;
}

const DataReflectionConditions& putInformationOnReflectionConditions(const BravaisType& brav_type, const Int4& irc_type)
{
//    static const DataReflectionConditions DATA_NONE("None", "", &is_not_extinct_none);
    static const DataReflectionConditions DATA_CUBIC_F[DATA_NUM_CUBIC_F]
	= {
			DataReflectionConditions("No condition:196,202,209,216,225" , "", &is_not_extinct_none),
			DataReflectionConditions("210","h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_4h00),
			DataReflectionConditions("203,227","0kl:k+l=4n,h0l:h+l=4n,hk0:h+k=4n,h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_4h00_40kl),
			DataReflectionConditions("219,226","hhl:h,l=2n,hkh:h,k=2n,hkk:h,k=2n", &is_not_extinct_2hhl),
			DataReflectionConditions("228","0kl:k+l=4n,h0l:h+l=4n,hk0:h+k=4n,hhl:h,l=2n,hkh:h,k=2n,hkk:h,k=2n,h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_40kl_2hhl_4h00),
			// A, face
			DataReflectionConditions("203f(x,0,0),203b(1/2,1/2,1/2),203a(0,0,0),210f(x,0,0),210b(1/2,1/2,1/2),210a(0,0,0),227f(x,0,0),227b(1/2,1/2,1/2),227a(0,0,0)"
			, "hkl:h=2n+1 or h+k+l=4n", &special_reflection_conditions_3h_4phkl),
			// B, face
			DataReflectionConditions("203d(5/8,5/8,5/8),203c(1/8,1/8,1/8),210d(5/8,5/8,5/8),210c(1/8,1/8,1/8),227d(5/8,5/8,5/8),227c(1/8,1/8,1/8)"
			, "hkl:h=2n+1 or h,k,l=4n+2 or h,k,l=4n", &special_reflection_conditions_3h_6ahkl_4ahkl),
    };

    static const DataReflectionConditions DATA_CUBIC_I[DATA_NUM_CUBIC_I]
	= {
			DataReflectionConditions("No condition:197,199,204,211,217,229" , "", &is_not_extinct_none),
			DataReflectionConditions( "206","0kl:k,l=2n,h0l:h,l=2n,hk0:h,k=2n", &is_not_extinct_20kl ),
			DataReflectionConditions( "214","h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_4h00 ),
			DataReflectionConditions( "220,230","hhl:2h+l=4n,hkh:2h+k=4n,hkk:h+2k=4n,h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_4hhl_4h00 ),
			//F
			DataReflectionConditions("220c(x,x,x),230e(x,x,x)" , "hkl:h=2n+1 or k=2n+1 or l=2n+1 or h+k+l=4n", &special_reflection_conditions_oddh_oddk_oddl_4phkl),
			//I
			DataReflectionConditions("214b(7/8,7/8,7/8),214a(1/8,1/8,1/8)" , "hkl:h=2n+1 or h,k,l=4n+2 or h,k,l=4n", &special_reflection_conditions_oddh_6ahkl_4ahkl),
			//J
			DataReflectionConditions("214d(5/8,0/1/4),214c(1/8,0,1/4)"
			, "hkl:h,k,l=2n,h+k+l=4n or h=8n,k=8n+4 and h+k+l=4n+2 or h,k=2n+1,l=4n+2 or h,k=8n+1,l=4n or h,k=8n+3,l=4n or h=8n+1,k=8n-1,l=4n or h=8n+3,k=8n-3,l=4n"
            , &special_reflection_conditions_2ahkl_4phkl_8h_12k_6phkl_oddahk_9ahk_4l_11ahk_4l_9h_7k_4l_11h),
			//J2
			DataReflectionConditions("220b(7/8,0,1/4),220a(3/8,0,1/4)"
			, "hkl:h,k=2n,h+k+l=4n or h=8n,k=8n+4 and h+k+l=4n+2 or h,k=2n+1,l=4n+2 or h=8n+1,k=8n+3,l=4n or h=8n+7,k=8n+5,l=4n or h=8n+1,k=8n+5,l=4n or h=8n+7,k=8n+3,l=4n"
            , &special_reflection_conditions_2hk_4phkl_8h_12k_6phkl_oddahk_6l_9h_11k_4l_15h_13k_9h_14k_4l_15h_11k_4l),
			//K
			DataReflectionConditions("214f(x,0,1/4)" , "hkl:h=2n+1 or h=4n,hhl:h=2n+1 or h+k+l=4n", &special_reflection_conditions_oddh_4h_oddh_4phkl_hhl),

			DataReflectionConditions("220d(x,0,1/4),230g(1/8,y,-y+1/4)" , "hkl:h=2n+1 or h=4n", &special_reflection_conditions_oddh_oddk_oddl_4hkl),
			//L
			DataReflectionConditions("230f(x,0,1/4)" , "hkl:2h+l=4n or 2k+h=4n or 2l+k=4n o r2l+h=4n or 2h+k=4n or 2k+l=4n"
		    , &special_reflection_conditions_4phl_4pkh_4plk_4plh_4phk_4pkl),
			//M
			DataReflectionConditions("230d(3/8,0,1/4),230c(1/8,0,1/4)" , "hkl:h,k=2n,h+k+l=4n or h,k=2n+1,l=4n+2 or h=8n,k=8n+4 and h+k+l=4n+2"
		    , &special_reflection_conditions_2ahk_4phkl_oddahk_6l_8h_12k_6phkl),
			//N
			DataReflectionConditions("230b(1/8,1/8,1/8)" , "hkl:h,k=2n+1,l=4n+2 or h,k,l=4n"
		    , &special_reflection_conditions_oddahk_6l_4ahkl),

    };
    static const DataReflectionConditions DATA_CUBIC_P[DATA_NUM_CUBIC_P]
	= {
			DataReflectionConditions("No condition:195,200,207,215,221" , "", &is_not_extinct_none),
			DataReflectionConditions( "198,208","h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_2h00 ),
			DataReflectionConditions( "201,224","h00:h=2n,0k0:k=2n,00l:l=2n,0kl:k+l=2n,h0l:h+l=2n,hk0:h+k=2n", &is_not_extinct_2h00_20kl ),
			DataReflectionConditions( "205","hk0:h=2n,0kl:k=2n,h0l:l=2n,h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_2hk0_2h00 ),
			DataReflectionConditions( "205(mirror-reversed)","hk0:k=2n,0kl:l=2n,h0l:h=2n,h00:h=2n,0k0:k=2n,00l:l=2n",&is_not_extinct_2hk0mirror_2h00),
			DataReflectionConditions( "212,213","h00:h=4n,0k0:k=4n,00l:l=4n", &is_not_extinct_4h00 ),
			DataReflectionConditions( "218,223","hhl:l=2n,hkh:k=2n,hkk:h=2n,h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_2hhl_2h00 ),
			DataReflectionConditions( "222","0kl:k+l=2n,h0l:h+l=2n,hk0:h+k=2n,hhl:l=2n,hkh:k=2n,hkk:h=2n,h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_2hhl_2h00_20kl ),
			//F
			DataReflectionConditions("208j(x,1/2,0),208i(x,0,1/2)"
			, "hkl:h=2n or k=2n or l=2n,hhl:l=2n,hkh:k=2n,hkk:h=2n", &special_reflection_conditions_2h_2k_2hhl_2hkh_2hkk),
			DataReflectionConditions("218h(x,0,1/2),218g(x,1/2,0),223j(1/4,y,y+1/2),223h(x,1/2,0),223g(x,0,1/2)"
			, "hkl:h=2n or k=2n or l=2n", &special_reflection_conditions_2h_2k_2lhkl),
			//G
			DataReflectionConditions("208f(1/4,1/2,0),208e(1/4,0,1/2),218d(1/4,0,1/2),218c(1/4,1/2,0),223d(1/4,1/2,0),223c(1/4,0,1/2)"
			, "hkl:h+k+l=2n or either one of h,k,l is 2n+1,4n and 4n+2", &special_reflection_conditions_2phkl_oddhkl_4hkl_6hkl),
			//H
			DataReflectionConditions("212b(5/8,5/8,5/8),212a(1/8,1/8,1/8),213b(7/8,7/8,7/8),213a(3/8,3/8,3/8)"
			, "hkl:two of h,k,l are odd or either one of h,k,l is 2n+1,k=4n and l=4n+2 or h,k,l=4n+2 or h,k,l=4n", &special_reflection_conditions_22hkl_2oddhkl_4k_6l_6ahkl_4hkl),



    };
    static const DataReflectionConditions DATA_HEXAGONAL[DATA_NUM_HEXAGONAL]
	= {
	        //Hexagonal :hex
			DataReflectionConditions("No condition:143,147,149,150,156,157,162,164,168,174,175,177,183,187,189,191" , "", &is_not_extinct_none),
			DataReflectionConditions( "144,145,151,152,153,154,171,172,180,181","000l:l=3n", &is_not_extinct_tr_3000l ),
			DataReflectionConditions( "158,165,185,188,193","h-h0l:l=2n,h0-hl:l=2n,0h-hl:l=2n", &is_not_extinct_tr_2hmh0l ),
			DataReflectionConditions( "159,163,186,190,194","hh-2hl:l=2n,h-2hhl:l=2n,-2hhhl:l=2n", &is_not_extinct_2hhm2hl ),
			DataReflectionConditions( "169,170,178,179","000l:l=6n", &is_not_extinct_hex_6000l ),
			DataReflectionConditions( "173,176,182","000l:l=2n", &is_not_extinct_hex_2000l ),
			DataReflectionConditions( "184,192","hh-2hl:l=2n,h-2hhl:l=2n,-2hhhl:l=2n,h-h0l:l=2n,h0-hl:l=2n,0h-hl:l=2n", &is_not_extinct_hex_2hhm2hl_2hmh0l_2000l ),
			//D
			DataReflectionConditions("159b(1/3,2/3,z),163f(1/3,2/3,z),163d(2/3,1/3,1/4),163c(1/3,2/3,1/4),173b(1/3,2/3,z),176f(1/3,2/3,z),176d(2/3,1/3,1/4),176c(1/3,2/3,1/4),182f(1/3,2/3,z),182d(1/3,2/3,3/4),182c(1/3,2/3,1/4),186b(1/3,2/3,z),190f(1/3,2/3,z),190d(2/3,1/3,1/4),190c(1/3,2/3,1/4),194f(1/3,2/3,z),194d(1/3,2/3,3/4),194c(1/3,2/3,1/4)"
			, "hkil:l=2n or h-k=3n+1 or h-k=3n+2", &special_reflection_conditions_2l_3shk_3s2hk),
			//E
			DataReflectionConditions("171b(1/2,1/2,z),172b(1/2,1/2,z),180f(1/2,0,z),180d(1/2,0,1/2),180c(1/2,0,0),181f(1/2,0,z),181d(1/2,0,1/2),181c(1/2,0,0)"
			, "hkil:h=2n+1 or k=2n+1 or l=3n", &special_reflection_conditions_oddh_oddk_3lhkil)
    };
    static const DataReflectionConditions DATA_TETRAGONAL_P[DATA_NUM_TETRAGONAL_P]
	= {
			DataReflectionConditions("No condition:75,81,83,89,99,111,115,123" , "", &is_not_extinct_none),
			DataReflectionConditions( "76,78,91,95","00l:l=4n", &is_not_extinct_400l ),
			DataReflectionConditions( "77,84,93","00l:l=2n", &is_not_extinct_200l ),
			DataReflectionConditions( "86","hk0:h+k=2n,00l:l=2n", &is_not_extinct_2hk0_200l ),
			DataReflectionConditions( "90,113","h00:h=2n,0k0:k=2n", &is_not_extinct_2h00_20k0 ),
			DataReflectionConditions( "92,96","00l:l=4n,h00:h=2n,0k0:k=2n", &is_not_extinct_400l_2h00_0k0 ),
			DataReflectionConditions( "94","h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_2h00 ),
			DataReflectionConditions( "100,117,127","h0l:l=2n,0kl:l=2n,h00:h=2n,0k0:k=2n", &is_not_extinct_t_2h0l_2h00 ),//t:Tetragonal
			DataReflectionConditions( "101,116,132","h0l:l=2n,0kl:l=2n", &is_not_extinct_t_2h0l ),
			DataReflectionConditions( "102,118,136","h0l:h+l=2n,0kl:k+l=2n,h00:h=2n,0k0:k=2n", &is_not_extinct_t_2ph0l_2h00 ),//p:plush+l=0
			DataReflectionConditions( "103,124","h0l:l=2n,0kl:l=2n,hhl:l=2n", &is_not_extinct_t_2h0l_2hhl ),
			DataReflectionConditions( "105,112,131","hhl:l=2n", &is_not_extinct_t_2hhl ),
			DataReflectionConditions( "106,135","h0l:h=2n,0kl:k=2n,hhl:l=2n,h00:h=2n,0k0:k=2n", &is_not_extinct_t_2h00_2h0l_2hhl ),
			DataReflectionConditions( "114","hhl:l=2n,h00:h=2n,0k0:k=2n", &is_not_extinct_t_2h00_2hhl ),
			DataReflectionConditions( "125","hk0:h+k=2n,h0l:h=2n,0kl:k=2n", &is_not_extinct_t_2phk0_2h0l ),
			DataReflectionConditions( "126","hk0:h+k=2n,h0l:h+l=2n,0kl:k+l=2n,hhl:l=2n", &is_not_extinct_t_2phk0_2hhl ),
			DataReflectionConditions( "104,128","h0l:h+l=2n,0kl:k+l=2n,hhl:l=2n,h00:h=2n,0k0:k=2n", &is_not_extinct_t_2ph0l_2hhl_2h00 ),
			DataReflectionConditions( "85,129","hk0:h+k=2n", &is_not_extinct_t_2phk0 ),
			DataReflectionConditions( "130","hk0:h+k=2n,h0l:l=2n,0kl:l=2n,hhl:l=2n", &is_not_extinct_t_2phk0_2h0l_2hll ),
			DataReflectionConditions( "133","hk0:h+k=2n,hhl:l=2n", &is_not_extinct_t_2phk0_2hll ),
			DataReflectionConditions( "134","hk0:h+k=2n,h0l:h+l=2n,0kl:k+l=2n", &is_not_extinct_t_23phk0 ),//3:hk0,0kl,h0l
			DataReflectionConditions( "137","hk0:h+k=2n,hhl:l=2n", &is_not_extinct_t_21phk0_2hhl ),//1:hk0,no h0k or 0kl
			DataReflectionConditions( "138","hk0:h+k=2n,h0l:l=2n,0kl:l=2n", &is_not_extinct_t_21phk0_2h0l ),//if Tetragonal and 2 conditions (h0l,0kl)no 1 or 3
	};
    static const DataReflectionConditions DATA_TETRAGONAL_I[DATA_NUM_TETRAGONAL_I]
	= {
			DataReflectionConditions("No condition:79,82,87,97,107,119,121,139" , "", &is_not_extinct_none),
			DataReflectionConditions( "80,98","00l:l=4n", &is_not_extinct_400l ),
			DataReflectionConditions( "88","hk0:h,k=2n,00l:l=4n", &is_not_extinct_t_21hk0_200l ),
			DataReflectionConditions( "108,120,140","h0l:h,l=2n,0kl:k,l=2n", &is_not_extinct_t_2ah0l ),//a:and h,l=2n
			DataReflectionConditions( "109,122","hhl:2h+l=4n,h-h0:h=2n", &is_not_extinct_t_4phhl_2hmh0 ),//mh:h-h0
			DataReflectionConditions( "110","hhl:2h+l=4n,h-h0:h=2n,h0l:h,l=2n,0kl:k,l=2n", &is_not_extinct_t_4phhl_2hmh0_2ah0l ),
			DataReflectionConditions( "141","hk0:h,k=2n,hhl:2h+l=4n", &is_not_extinct_t_2hk0_4phhl ),
			DataReflectionConditions( "142","hk0:h,k=2n,hhl:2h+l=4n,h0l:h,l=2n,0kl:k,l=2n", &is_not_extinct_t_23ahk0_4phhl ),

			//A,body
			DataReflectionConditions("80a(0,0,z),88e(0,0,z),88b(0,0,1/2),88a(0,0,0),98c(0,0,z),98b(0,0,1/2),98a(0,0,0),109a(0,0,z),122c(0,0,z),122b(0,0,1/2),122a(0,0,0),141g(x,x,0),141e(0,0,z),141b(0,0,1/2),141a(0,0,0),142f(x,x,1/4)"
			, "hkl:h=2n+1 or 2h+l=4n", &special_reflection_conditions_3l_4p2hl),
			//B,body
			DataReflectionConditions("88d(0,1/4,5/8),88c(0,1/4,1/8),141d(0,1/4,5/8),141c(0,1/4,1/8)"
			, "hkl:l=2n+1 or h,k=2n,h+k+l=4n", &special_reflection_conditions_3l_2ahk_4phkl),
			//C
			DataReflectionConditions("141f(x,1/4,1/8),142e(1/4,y,1/8)"
			, "hkl:l=2n+1 or h=2n", &special_reflection_conditions_3l_2h),
	};
    static const DataReflectionConditions DATA_RHOMBOHEDRAL_RHOM_AXIS[DATA_NUM_RHOMBOHEDRAL_RHOM_AXIS]
    = {
    		DataReflectionConditions("No condition:146,148,155,160,166" , "", &is_not_extinct_none),
    		DataReflectionConditions("161,167","hhl:l=2n,hkh:k=2n,hkk:h=2n", &is_not_extinct_rhom_2hhl ),
    };

    static const DataReflectionConditions DATA_RHOMBOHEDRAL_HEX_AXIS[DATA_NUM_RHOMBOHEDRAL_HEX_AXIS]
    = {
    		DataReflectionConditions("No condition:146,148,155,160,166" , "", &is_not_extinct_none),
    		DataReflectionConditions("161,167", "h-h0l:l=2n,h0-hl:l=2n,0h-hl:l=2n", &is_not_extinct_rhom_hex_2hmh0l),
    };

    //Monoclinic P
    static const DataReflectionConditions DATA_MONOCLINIC_P_A_AXIS[DATA_NUM_MONOCLINIC_P_A_AXIS]
    = {
    		DataReflectionConditions("No condition:3,6,10" , "", &is_not_extinct_none),
    		DataReflectionConditions( "4,11","h00:h=2n", &is_not_extinct_mono_2h00 ),
    		DataReflectionConditions( "7,13","0kl:k=2n", &is_not_extinct_mono_20kl ),
    		DataReflectionConditions( "7(cell choice 2),13(cell choice 2)","0kl:k+l=2n", &is_not_extinct_mono_2p0kl ),
    		DataReflectionConditions( "7(cell choice 3),13(cell choice 3)","0kl:l=2n", &is_not_extinct_mono_2l0kl ),//l in 2l:l=2n
    		DataReflectionConditions( "14","0kl:k=2n,h00:h=2n", &is_not_extinct_mono_20kl_2h00 ),
    		DataReflectionConditions( "14(cell choice 2)","0kl:k+l=2n,0k0:k=2n", &is_not_extinct_mono_2p0kl_20k0 ),
    		DataReflectionConditions( "14(cell choice 3)","0kl:l=2n,h00:h=2n", &is_not_extinct_mono_2d0kl_2h00 )
     };


    static const DataReflectionConditions DATA_MONOCLINIC_P_B_AXIS[DATA_NUM_MONOCLINIC_P_B_AXIS]
    ={
    		DataReflectionConditions("No condition:3,6,10" , "", &is_not_extinct_none),
    		DataReflectionConditions( "4,11","0k0:k=2n", &is_not_extinct_mono_20k0 ),
    		DataReflectionConditions( "7,13","h0l:l=2n", &is_not_extinct_mono_2h0l ),
    		DataReflectionConditions( "7(cell choice 2),13(cell choice 2)","h0l:h+l=2n", &is_not_extinct_mono_2ph0l ),
			DataReflectionConditions( "7(cell choice 3),13(cell choice 3)","h0l:h=2n", &is_not_extinct_mono_2hh0l ),
			DataReflectionConditions( "14","h0l:l=2n,0k0:h=2n", &is_not_extinct_mono_2h0l_20k0 ),
			DataReflectionConditions( "14(cell choice 2)","h0l:h+l=2n,0k0:k=2n", &is_not_extinct_mono_2ph0l_20k0 ),
			DataReflectionConditions( "14(cell choice 3)","h0l:h=2n,0k0:l=2n", &is_not_extinct_mono_2dh0l_20k0 )
     };

    static const DataReflectionConditions DATA_MONOCLINIC_P_C_AXIS[DATA_NUM_MONOCLINIC_P_C_AXIS]
    ={
    		DataReflectionConditions("No condition:3,6,10" , "", &is_not_extinct_none),
    		DataReflectionConditions( "4,11","00l:l=2n", &is_not_extinct_mono_200l ),
    		DataReflectionConditions( "7,13","hk0:h=2n", &is_not_extinct_mono_2hk0 ),
			DataReflectionConditions( "7(cell choice 2),13(cell choice 2)","hk0:h+k=2n", &is_not_extinct_mono_2phk0 ),
			DataReflectionConditions( "7(cell choice 3),13(cell choice 3)","hk0:k=2n", &is_not_extinct_mono_2khk0 ),
			DataReflectionConditions( "14","hk0:h=2n,00l:l=2n", &is_not_extinct_mono_2hk0_200l ),
			DataReflectionConditions( "14(cell choice 2)","hk0:h+k=2n,00l:k=2n", &is_not_extinct_mono_2phk0_200l ),
			DataReflectionConditions( "14(cell choice 3)","hk0:k=2n,00l:l=2n", &is_not_extinct_mono_2dhk0_200l )
    };

    static const DataReflectionConditions DATA_MONOCLINIC_B_A_AXIS[DATA_NUM_MONOCLINIC_B_A_AXIS]
    = {
    		DataReflectionConditions("No condition:5,8,12" , "", &is_not_extinct_none),
    		DataReflectionConditions( "9,15", "0kl:k=2n", &is_not_extinct_mono_20kl )
    };

     static const DataReflectionConditions DATA_MONOCLINIC_B_B_AXIS[DATA_NUM_MONOCLINIC_B_B_AXIS]
    = {
    		DataReflectionConditions("No condition:5,8,12" , "", &is_not_extinct_none),
        	DataReflectionConditions( "9,15", "h0l:l=2n", &is_not_extinct_mono_2h0l )
    };

    static const DataReflectionConditions DATA_MONOCLINIC_B_C_AXIS[DATA_NUM_MONOCLINIC_B_C_AXIS]
    = {
    		DataReflectionConditions("No condition:5,8,12" , "", &is_not_extinct_none),
        	DataReflectionConditions( "9,15", "hk0:h=2n", &is_not_extinct_mono_2hk0 )
    };
    static const DataReflectionConditions DATA_ORTHORHOMBIC_P[DATA_NUM_ORTHORHOMBIC_P]
    = {
    		DataReflectionConditions("No condition:16,25,47" , "", &is_not_extinct_none),
			DataReflectionConditions("19","h00:h=2n,0k0:k=2n,00l:l=2n", &is_not_extinct_orth_2dh00),
			DataReflectionConditions("48","hk0:h+k=2n,h0l:h+l=2n,0kl:k+l=2n", &is_not_extinct_orth_23phk0, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("61","0kl:k=2n,h0l:l=2n,hk0:h=2n", &is_not_extinct_orth_2d30kl),//3:0kl,h0l,hk0
			DataReflectionConditions("61(mirror-reversed)","0kl:l=2n,h0l:h=2n,hk0:k=2n", &is_not_extinct_orth_2d30kl, DataReflectionConditions::AXIS_ACB),//3:0kl,h0l,hk0

			DataReflectionConditions("17","00l:l=2n", &standard_function_for_abc_200l, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("17(cab)","h00:h=2n", &standard_function_for_abc_200l, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("17(bca)","0k0:k=2n", &standard_function_for_abc_200l, DataReflectionConditions::AXIS_BCA),

			DataReflectionConditions("18","h00:h=2n,0k0:k=2n", &standard_function_for_abc_22h00, DataReflectionConditions::AXIS_ABC),//2:h00,0k0
			DataReflectionConditions("18(cab)","0k0:k=2n,00l:l=2n", &standard_function_for_abc_22h00, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("18(bca)","00l:l=2n,h00:h=2n", &standard_function_for_abc_22h00, DataReflectionConditions::AXIS_BCA),

			DataReflectionConditions("27,49","0kl:l=2n,h0l:l=2n", &standard_function_for_abc_220kl, DataReflectionConditions::AXIS_ABC),//2:0kl,h0l
			DataReflectionConditions("27(cab),49(cab)","h0l:h=2n,hk0:h=2n", &standard_function_for_abc_220kl, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("27(bca),49(bca)","hk0:k=2n,0kl:k=2n", &standard_function_for_abc_220kl, DataReflectionConditions::AXIS_BCA),

			DataReflectionConditions("32,55","0kl:k=2n,h0l:h=2n", &standard_function_for_abc_2d20kl, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("32(cab),55(cab)","h0l:l=2n,hk0:k=2n", &standard_function_for_abc_2d20kl, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("32(bca),55(bca)","hk0:h=2n,0kl:l=2n", &standard_function_for_abc_2d20kl, DataReflectionConditions::AXIS_BCA),

			DataReflectionConditions("34,58","0kl:k+l=2n,h0l:h+l=2n", &standard_function_for_abc_2p20kl, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("34(cab),58(cab)","h0l:h+l=2n,hk0:h+k=2n", &standard_function_for_abc_2p20kl, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("34(bca),58(bca)","hk0:k+h=2n,0kl:k+l=2n", &standard_function_for_abc_2p20kl, DataReflectionConditions::AXIS_ABC),

			DataReflectionConditions("50","0kl:k=2n,h0l:h=2n,hk0:h+k=2n", &standard_function_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("50(cab)","h0l:l=2n,hk0:k=2n,0kl:k+l=2n", &standard_function_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("50(bca)","hk0:h=2n,0kl:l=2n,h0l:h+l=2n", &standard_function_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_BAC),


			DataReflectionConditions("56","0kl:l=2n,h0l:l=2n,hk0:h+k=2n", &standard_function_for_abc, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("56(cab)","h0l:h=2n,hk0:h=2n,0kl:k+l=2n", &standard_function_for_abc, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("56(bca)","hk0:k=2n,0kl:k=2n,h0l:h+l=2n", &standard_function_for_abc, DataReflectionConditions::AXIS_BCA),

			// case of (abc), the string conditions are different, the function are the same
            DataReflectionConditions("26,28(-cba),51(bca)",      "h0l:l=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_ABC),
            DataReflectionConditions("26(cab),28(a-cb),51",      "hk0:h=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_CAB),
            DataReflectionConditions("26(bca),28(ba-c),51(cab)", "0kl:k=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_BCA),
            DataReflectionConditions("26(ba-c),28(bca),51(-cba)", "0kl:l=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_BAC),
            DataReflectionConditions("26(-cba),28,51(a-cb)",      "h0l:h=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_CBA),
            DataReflectionConditions("26(a-cb),28(cab),51(ba-c)", "hk0:k=2n", &standard_function2_for_abc, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("29,57" ,           "0kl:l=2n,h0l:h=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("29(cab),57(cab)" , "h0l:h=2n,hk0:k=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("29(bca),57(bca)" , "h0l:k=2n,0kl:l=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("29(ba-c),57(ba-c)" , "h0l:l=2n,0kl:k=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("29(-cba),57(-cba)" , "hk0:h=2n,h0l:l=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("29(a-cb),57(a-cb)" , "0kl:k=2n,hk0:h=2n", &standard_function2_for_abc_22d0kl, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("30,53" ,           "0kl:k+l=2n,h0l:l=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("30(cab),53(cab)" , "h0l:h+l=2n,hk0:h=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("30(bca),53(bca)" , "hk0:h+k=2n,0kl:k=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("30(ba-c),53(ba-c)" , "h0l:h+l=2n,0kl:l=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("30(-cba),53(-cba)" , "hk0:h+k=2n,h0l:h=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("30(a-cb),53(a-cb)" , "0kl:k+l=2n,hk0:k=2n", &standard_function2_for_abc_2p0kl_2h0l, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("31,59(bca)" , "h0l:h+l=2n", &standard_function2_for_abc_2ph0l, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("31(cab),59" , "hk0:h+k=2n", &standard_function2_for_abc_2ph0l, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("31(bca),59(cab)", "0kl:k+l=2n", &standard_function2_for_abc_2ph0l, DataReflectionConditions::AXIS_BCA),

			DataReflectionConditions("33,62" ,           "0kl:k+l=2n,h0l:h=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("33(cab),62(cab)" , "h0l:h+l=2n,hk0:k=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("33(bca),62(bca)" , "hk0:h+k=2n,0kl:l=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("33(ba-c),62(ba-c)" , "h0l:h+l=2n,0kl:k=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("33(-cba),62(-cba)" , "hk0:h+k=2n,h0l:l=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("33(a-cb),62(a-cb)" , "0kl:k+l=2n,hk0:h=2n", &standard_function2_for_abc_2p0kl_2dh0l, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("52" ,      "0kl:k+l=2n,h0l:h+l=2n,hk0:h=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("52(cab)" , "h0l:h+l=2n,hk0:h+k=2n,0kl:k=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("52(bca)" , "hk0:h+k=2n,0kl:k+l=2n,h0l:l=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("52(ba-c)" , "h0l:h+l=2n,0kl:k+l=2n,hk0:k=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("52(-cba)" , "hk0:h+k=2n,h0l:h+l=2n,0kl:l=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("52(a-cb)" , "0kl:k+l=2n,hk0:h+k=2n,h0l:h=2n", &standard_function2_for_abc_22p0kl_2hk0, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("54" ,      "0kl:l=2n,h0l:l=2n,hk0:h=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("54(cab)" , "h0l:h=2n,hk0:h=2n,0kl:k=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("54(bca)" , "hk0:k=2n,0kl:k=2n,h0l:l=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("54(ba-c)" , "h0l:l=2n,0kl:l=2n,hk0:k=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("54(-cba)" , "hk0:h=2n,h0l:h=2n,0kl:l=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("54(a-cb)" , "0kl:k=2n,hk0:k=2n,h0l:h=2n", &standard_function2_for_abc_230kl, DataReflectionConditions::AXIS_ACB),

			DataReflectionConditions("60" ,      "0kl:k=2n,h0l:l=2n,hk0:h+k=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_ABC),
			DataReflectionConditions("60(cab)" , "h0l:l=2n,hk0:h=2n,0kl:k+l=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_CAB),
			DataReflectionConditions("60(bca)" , "hk0:h=2n,0kl:k=2n,h0l:h+l=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_BCA),
			DataReflectionConditions("60(ba-c)" , "h0l:h=2n,0kl:l=2n,hk0:h+k=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_BAC),
			DataReflectionConditions("60(-cba)" , "hk0:k=2n,h0l:h=2n,0kl:k+l=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_CBA),
			DataReflectionConditions("60(a-cb)" , "0kl:l=2n,hk0:k=2n,h0l:h+l=2n", &standard_function2_for_abc_220kl_2phk0, DataReflectionConditions::AXIS_ACB),
    };
    static const DataReflectionConditions DATA_ORTHORHOMBIC_C[DATA_NUM_ORTHORHOMBIC_C]
    = {
    		DataReflectionConditions("No condition:21,35,38,65" , "", &is_not_extinct_none),
    		DataReflectionConditions("20" , "00l:l=2n", &standard_function_for_abc_200l),
    		DataReflectionConditions("36,63" , "h0l:l=2n", &standard_function2_for_abc),
			DataReflectionConditions("40(bca)" , "0kl:l=2n", &standard_function_for_abc_210kl),
    		DataReflectionConditions("37,66" , "0kl:l=2n,h0l:l=2n", &standard_function_for_abc_220kl),
    		DataReflectionConditions("39(bca),67" , "hk0:h,k=2n", &standard_function_for_abc_2ahk0),
    		DataReflectionConditions("41(bca),64(ba-c)" , "hk0:h,k=2n,0kl:l=2n", &standard_function_for_abc_2ahk0_20kl, DataReflectionConditions::AXIS_BAC),
    		DataReflectionConditions("41(-cba),64" , "hk0:h,k=2n,h0l:l=2n", &standard_function_for_abc_2ahk0_2h0l, DataReflectionConditions::AXIS_ABC),
    		DataReflectionConditions("68" , "0kl:l=2n,h0l:l=2n,hk0:h,k=2n", &standard_function_for_abc_220kl_2ahk0),

    };
    static const DataReflectionConditions DATA_ORTHORHOMBIC_F[DATA_NUM_ORTHORHOMBIC_F]
    = {
    		DataReflectionConditions("No condition:22,42,69" , "", &is_not_extinct_none),
    		DataReflectionConditions("43" , "0kl:k+l=4n,h0l:h+l=4n", &standard_function_for_abc_42p0kl, DataReflectionConditions::AXIS_ABC),
    		DataReflectionConditions("43(cab)" , "h0l:h+l=4n,hk0:h+k=4n", &standard_function_for_abc_42p0kl, DataReflectionConditions::AXIS_CAB),
    		DataReflectionConditions("43(bca)" , "hk0:h+k=4n,0kl:k+l=4n", &standard_function_for_abc_42p0kl, DataReflectionConditions::AXIS_BCA),

    		DataReflectionConditions("70" , "0kl:k+l=4n,h0l:h+l=4n,hk0:h+k=4n", &standard_function_for_abc_43p0kl),

			// A
			DataReflectionConditions("43a(0,0,z),70g(0,0,z),70f(0,y,0),70e(x,0,0),70b(0,0,1/2),70a(0,0,0)" , "hkl:h=2n+1 or h+k+l=4n", &special_reflection_conditions_3h_4phkl),

			// B
			DataReflectionConditions("70d(5/8,5/8,5/8),70c(1/8,1/8,1/8)" , "hkl:h=2n+1 or h,k,l=4n+2 or h,k,l=4n", &special_reflection_conditions_3h_6ahkl_4ahkl),
    };

    static const DataReflectionConditions DATA_ORTHORHOMBIC_I[DATA_NUM_ORTHORHOMBIC_I]
    = {
    		DataReflectionConditions("No condition:23,24,44,71" , "", &is_not_extinct_none),
    		DataReflectionConditions("45,72" , "0kl:k,l=2n,h0l:h,l=2n", &standard_function_for_abc_2a0kl, DataReflectionConditions::AXIS_ABC),
    		DataReflectionConditions("45(cab),72(cab)" , "h0l:h,l=2n,hk0:h,k=2n", &standard_function_for_abc_2a0kl, DataReflectionConditions::AXIS_CAB),
    		DataReflectionConditions("45(bca),72(bca)" , "hk0:h,k=2n,0kl:k,l=2n", &standard_function_for_abc_2a0kl, DataReflectionConditions::AXIS_BCA),

    		DataReflectionConditions("46,74" , "h0l:h,l=2n", &standard_function_for_abc_2ah0l, DataReflectionConditions::AXIS_ABC),
    		DataReflectionConditions("46(cab),74(cab)" , "hk0:h,k=2n", &standard_function_for_abc_2ah0l, DataReflectionConditions::AXIS_CAB),
    		DataReflectionConditions("46(bca),74(bca)" , "0kl:k,l=2n", &standard_function_for_abc_2ah0l, DataReflectionConditions::AXIS_BCA),

    		DataReflectionConditions("73" , "0kl:k,l=2n,h0l:h,l=2n,hk0:h,k=2n", &standard_function_for_abc_23a0kl),

    };
    static const DataReflectionConditions DATA_TRICLINIC[DATA_NUM_TRICLINIC]
    = {
    		DataReflectionConditions("No condition:1,2" , "", &is_not_extinct_none),
    };

//    if( irc_type < 0 ){ return DATA_NONE; }

    if( brav_type.enumBravaisType() == Cubic_F )
    {
        return DATA_CUBIC_F[(size_t) irc_type];
    }
    if( brav_type.enumBravaisType() == Cubic_I )
    {
        return DATA_CUBIC_I[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Cubic_P )
    {
        return DATA_CUBIC_P[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Hexagonal )
    {
        return DATA_HEXAGONAL[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Tetragonal_P )
    {
        return DATA_TETRAGONAL_P[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Tetragonal_I )
    {
        return DATA_TETRAGONAL_I[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Rhombohedral )
    {
            if( brav_type.enumRHaxis() == Rho_Axis )
            {
                    return DATA_RHOMBOHEDRAL_RHOM_AXIS[(size_t) irc_type];
            }
            if( brav_type.enumRHaxis() == Hex_Axis )
            {
                    return DATA_RHOMBOHEDRAL_HEX_AXIS[(size_t) irc_type];
            }
     }
    if ( brav_type.enumBravaisType() == Orthorhombic_P )
    {
        return DATA_ORTHORHOMBIC_P[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Orthorhombic_C )
    {
        return DATA_ORTHORHOMBIC_C[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Orthorhombic_F )
    {
        return DATA_ORTHORHOMBIC_F[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Orthorhombic_I )
    {
        return DATA_ORTHORHOMBIC_I[(size_t) irc_type];
    }
    if ( brav_type.enumBravaisType() == Monoclinic_P )
    {
        if( brav_type.enumABCaxis() == A_Axis )
        {
                return DATA_MONOCLINIC_P_A_AXIS[(size_t) irc_type];
        }
        if( brav_type.enumABCaxis() == B_Axis )
        {
                return DATA_MONOCLINIC_P_B_AXIS[(size_t) irc_type];
        }
        if( brav_type.enumABCaxis() == C_Axis )
        {
                return DATA_MONOCLINIC_P_C_AXIS[(size_t) irc_type];
        }
    }
    if ( brav_type.enumBravaisType() == Monoclinic_B )
    {
        if( brav_type.enumABCaxis() == A_Axis )
        {
                return DATA_MONOCLINIC_B_A_AXIS[(size_t) irc_type];
        }
        if( brav_type.enumABCaxis() == B_Axis )
        {
                return DATA_MONOCLINIC_B_B_AXIS[(size_t) irc_type];
        }
        if( brav_type.enumABCaxis() == C_Axis )
        {
                return DATA_MONOCLINIC_B_C_AXIS[(size_t) irc_type];
        }
    }
    if ( brav_type.enumBravaisType() == Triclinic )
    {
        return DATA_TRICLINIC[(size_t) irc_type];
    }

    assert( false );
    return DATA_TRICLINIC[0];
}

string DataReflectionConditions::putShortStringType() const
{
	string ans;
	istringstream iss(type);
	ZErrorMessage zerr = getdelim(iss, ans, "(");
	if( zerr.putErrorType() == ZErrorDelimiterNotFound || isalpha(*(ans.rbegin()) ) ) return ans;
	return type;
}
