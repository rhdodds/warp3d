c     ****************************************************************
c     *                                                              *
c     *                      subroutine efrwrd_gen                   *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/27/09 rhd (threads)     *
c     *                                                              *
c     *     this subroutine performs forward reduction for a         *
c     *     block of non-conflicting, similar elements.              * 
c     *                                                              *
c     *     this version is the general case which handles any       *
c     *     element type.  However, it is not as fast as the         *
c     *     unrolled version, so is provided for debugging.          *
c     *                                                              *
c     ****************************************************************
c     
c     
c     
      subroutine efrwrd_gen( span, totdof, edest, z, pcm, utsz )
      implicit integer (a-z)
$add param_def
      dimension edest(totdof,*)
#sgl      real
#dbl      double precision
     &     pcm(utsz,*), z(*), ze(mxedof)
      dimension tops (mxedof)
c
      tops(1) = 1
      do i = 2, totdof
         tops(i) = tops(i-1) + i-1
      end do
c
c$OMP PARALLEL DO DEFAULT( shared )
c$OMP&            PRIVATE( i, j, k, ze )
      do i = 1, span
	 do j = 1, totdof
	    ze(j) = z(edest(j,i))
         end do
c
	 do j = 2, totdof
	    do k = 1, j - 1
	       ze(j) = ze(j) - pcm( tops(j)+k-1,i) * ze(k)
            end do
         end do
c
	 do j = 1, totdof
	    z(edest(j,i)) = ze(j) 
         end do
c
      end do
c$OMP END PARALLEL DO
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine efrwrd8                      *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/27/09 rhd (threads)     *
c     *                                                              *
c     *     this subroutine performs forward reduction for a         *
c     *     block of non-conflicting, similar elements.              * 
c     *                                                              *
c     *     this version is coded specifically for 8 node elements.  *
c     *                                                              *
c     ****************************************************************
c     
c     
c     
      subroutine efrwrd8( span, totdof, edest, z, pcm, utsz )
      implicit integer (a-z)
$add param_def
      dimension edest(totdof,*)
#sgl      real
#dbl      double precision
     &     pcm(utsz,*), z(*), 
     &     ze1,ze2,ze3,ze4,ze5,ze6,ze7,ze8,ze9,
     &     ze10,ze11,ze12,ze13,ze14,ze15,ze16,ze17,ze18,ze19,
     &     ze20,ze21,ze22,ze23,ze24
c
c$OMP PARALLEL DO PRIVATE( i,
c$OMP&         ze1,ze2,ze3,ze4,ze5,ze6,ze7,ze8,ze9,
c$OMP&         ze10,ze11,ze12,ze13,ze14,ze15,ze16,ze17,ze18,ze19,
c$OMP&         ze20,ze21,ze22,ze23,ze24 )
c
      do i = 1, span
         ze1 = z(edest(1,i))
         ze2 = z(edest(2,i))
         ze3 = z(edest(3,i))
         ze4 = z(edest(4,i))
         ze5 = z(edest(5,i))
         ze6 = z(edest(6,i))
         ze7 = z(edest(7,i))
         ze8 = z(edest(8,i))
         ze9 = z(edest(9,i))
         ze10 = z(edest(10,i))
         ze11 = z(edest(11,i))
         ze12 = z(edest(12,i))
         ze13 = z(edest(13,i))
         ze14 = z(edest(14,i))
         ze15 = z(edest(15,i))
         ze16 = z(edest(16,i))
         ze17 = z(edest(17,i))
         ze18 = z(edest(18,i))
         ze19 = z(edest(19,i))
         ze20 = z(edest(20,i))
         ze21 = z(edest(21,i))
         ze22 = z(edest(22,i))
         ze23 = z(edest(23,i))
         ze24 = z(edest(24,i))
c 
         ze2 = ze2 -
     &      pcm(   2,i)*ze1
c 
         ze3 = ze3 -
     &      pcm(   4,i)*ze1 - 
     &      pcm(   5,i)*ze2
c 
         ze4 = ze4 -
     &      pcm(   7,i)*ze1 - 
     &      pcm(   8,i)*ze2 - pcm(   9,i)*ze3
c 
         ze5 = ze5 -
     &      pcm(  11,i)*ze1 - 
     &      pcm(  12,i)*ze2 - pcm(  13,i)*ze3 - pcm(  14,i)*ze4
c 
         ze6 = ze6 -
     &      pcm(  16,i)*ze1 - 
     &      pcm(  17,i)*ze2 - pcm(  18,i)*ze3 - pcm(  19,i)*ze4 - 
     &      pcm(  20,i)*ze5
c 
         ze7 = ze7 -
     &      pcm(  22,i)*ze1 - 
     &      pcm(  23,i)*ze2 - pcm(  24,i)*ze3 - pcm(  25,i)*ze4 - 
     &      pcm(  26,i)*ze5 - pcm(  27,i)*ze6
c 
         ze8 = ze8 -
     &      pcm(  29,i)*ze1 - 
     &      pcm(  30,i)*ze2 - pcm(  31,i)*ze3 - pcm(  32,i)*ze4 - 
     &      pcm(  33,i)*ze5 - pcm(  34,i)*ze6 - pcm(  35,i)*ze7
c 
         ze9 = ze9 -
     &      pcm(  37,i)*ze1 - 
     &      pcm(  38,i)*ze2 - pcm(  39,i)*ze3 - pcm(  40,i)*ze4 - 
     &      pcm(  41,i)*ze5 - pcm(  42,i)*ze6 - pcm(  43,i)*ze7 - 
     &      pcm(  44,i)*ze8
c 
         ze10 = ze10 -
     &      pcm(  46,i)*ze1 - 
     &      pcm(  47,i)*ze2 - pcm(  48,i)*ze3 - pcm(  49,i)*ze4 - 
     &      pcm(  50,i)*ze5 - pcm(  51,i)*ze6 - pcm(  52,i)*ze7 - 
     &      pcm(  53,i)*ze8 - pcm(  54,i)*ze9
c 
         ze11 = ze11 -
     &      pcm(  56,i)*ze1 - 
     &      pcm(  57,i)*ze2 - pcm(  58,i)*ze3 - pcm(  59,i)*ze4 - 
     &      pcm(  60,i)*ze5 - pcm(  61,i)*ze6 - pcm(  62,i)*ze7 - 
     &      pcm(  63,i)*ze8 - pcm(  64,i)*ze9 - pcm(  65,i)*ze10
c 
         ze12 = ze12 -
     &      pcm(  67,i)*ze1 - 
     &      pcm(  68,i)*ze2 - pcm(  69,i)*ze3 - pcm(  70,i)*ze4 - 
     &      pcm(  71,i)*ze5 - pcm(  72,i)*ze6 - pcm(  73,i)*ze7 - 
     &      pcm(  74,i)*ze8 - pcm(  75,i)*ze9 - pcm(  76,i)*ze10 - 
     &      pcm(  77,i)*ze11
c 
         ze13 = ze13 -
     &      pcm(  79,i)*ze1 - 
     &      pcm(  80,i)*ze2 - pcm(  81,i)*ze3 - pcm(  82,i)*ze4 - 
     &      pcm(  83,i)*ze5 - pcm(  84,i)*ze6 - pcm(  85,i)*ze7 - 
     &      pcm(  86,i)*ze8 - pcm(  87,i)*ze9 - pcm(  88,i)*ze10 - 
     &      pcm(  89,i)*ze11 - pcm(  90,i)*ze12
c 
         ze14 = ze14 -
     &      pcm(  92,i)*ze1 - 
     &      pcm(  93,i)*ze2 - pcm(  94,i)*ze3 - pcm(  95,i)*ze4 - 
     &      pcm(  96,i)*ze5 - pcm(  97,i)*ze6 - pcm(  98,i)*ze7 - 
     &      pcm(  99,i)*ze8 - pcm( 100,i)*ze9 - pcm( 101,i)*ze10 - 
     &      pcm( 102,i)*ze11 - pcm( 103,i)*ze12 - pcm( 104,i)*ze13
c 
         ze15 = ze15 -
     &      pcm( 106,i)*ze1 - 
     &      pcm( 107,i)*ze2 - pcm( 108,i)*ze3 - pcm( 109,i)*ze4 - 
     &      pcm( 110,i)*ze5 - pcm( 111,i)*ze6 - pcm( 112,i)*ze7 - 
     &      pcm( 113,i)*ze8 - pcm( 114,i)*ze9 - pcm( 115,i)*ze10 - 
     &      pcm( 116,i)*ze11 - pcm( 117,i)*ze12 - pcm( 118,i)*ze13 - 
     &      pcm( 119,i)*ze14
c 
         ze16 = ze16 -
     &      pcm( 121,i)*ze1 - 
     &      pcm( 122,i)*ze2 - pcm( 123,i)*ze3 - pcm( 124,i)*ze4 - 
     &      pcm( 125,i)*ze5 - pcm( 126,i)*ze6 - pcm( 127,i)*ze7 - 
     &      pcm( 128,i)*ze8 - pcm( 129,i)*ze9 - pcm( 130,i)*ze10 - 
     &      pcm( 131,i)*ze11 - pcm( 132,i)*ze12 - pcm( 133,i)*ze13 - 
     &      pcm( 134,i)*ze14 - pcm( 135,i)*ze15
c 
         ze17 = ze17 -
     &      pcm( 137,i)*ze1 - 
     &      pcm( 138,i)*ze2 - pcm( 139,i)*ze3 - pcm( 140,i)*ze4 - 
     &      pcm( 141,i)*ze5 - pcm( 142,i)*ze6 - pcm( 143,i)*ze7 - 
     &      pcm( 144,i)*ze8 - pcm( 145,i)*ze9 - pcm( 146,i)*ze10 - 
     &      pcm( 147,i)*ze11 - pcm( 148,i)*ze12 - pcm( 149,i)*ze13 - 
     &      pcm( 150,i)*ze14 - pcm( 151,i)*ze15 - pcm( 152,i)*ze16
c 
         ze18 = ze18 -
     &      pcm( 154,i)*ze1 - 
     &      pcm( 155,i)*ze2 - pcm( 156,i)*ze3 - pcm( 157,i)*ze4 - 
     &      pcm( 158,i)*ze5 - pcm( 159,i)*ze6 - pcm( 160,i)*ze7 - 
     &      pcm( 161,i)*ze8 - pcm( 162,i)*ze9 - pcm( 163,i)*ze10 - 
     &      pcm( 164,i)*ze11 - pcm( 165,i)*ze12 - pcm( 166,i)*ze13 - 
     &      pcm( 167,i)*ze14 - pcm( 168,i)*ze15 - pcm( 169,i)*ze16 - 
     &      pcm( 170,i)*ze17
c 
         ze19 = ze19 -
     &      pcm( 172,i)*ze1 - 
     &      pcm( 173,i)*ze2 - pcm( 174,i)*ze3 - pcm( 175,i)*ze4 - 
     &      pcm( 176,i)*ze5 - pcm( 177,i)*ze6 - pcm( 178,i)*ze7 - 
     &      pcm( 179,i)*ze8 - pcm( 180,i)*ze9 - pcm( 181,i)*ze10 - 
     &      pcm( 182,i)*ze11 - pcm( 183,i)*ze12 - pcm( 184,i)*ze13 - 
     &      pcm( 185,i)*ze14 - pcm( 186,i)*ze15 - pcm( 187,i)*ze16 - 
     &      pcm( 188,i)*ze17 - pcm( 189,i)*ze18
c 
         ze20 = ze20 -
     &      pcm( 191,i)*ze1 - 
     &      pcm( 192,i)*ze2 - pcm( 193,i)*ze3 - pcm( 194,i)*ze4 - 
     &      pcm( 195,i)*ze5 - pcm( 196,i)*ze6 - pcm( 197,i)*ze7 - 
     &      pcm( 198,i)*ze8 - pcm( 199,i)*ze9 - pcm( 200,i)*ze10 - 
     &      pcm( 201,i)*ze11 - pcm( 202,i)*ze12 - pcm( 203,i)*ze13 - 
     &      pcm( 204,i)*ze14 - pcm( 205,i)*ze15 - pcm( 206,i)*ze16 - 
     &      pcm( 207,i)*ze17 - pcm( 208,i)*ze18 - pcm( 209,i)*ze19
c 
         ze21 = ze21 -
     &      pcm( 211,i)*ze1 - 
     &      pcm( 212,i)*ze2 - pcm( 213,i)*ze3 - pcm( 214,i)*ze4 - 
     &      pcm( 215,i)*ze5 - pcm( 216,i)*ze6 - pcm( 217,i)*ze7 - 
     &      pcm( 218,i)*ze8 - pcm( 219,i)*ze9 - pcm( 220,i)*ze10 - 
     &      pcm( 221,i)*ze11 - pcm( 222,i)*ze12 - pcm( 223,i)*ze13 - 
     &      pcm( 224,i)*ze14 - pcm( 225,i)*ze15 - pcm( 226,i)*ze16 - 
     &      pcm( 227,i)*ze17 - pcm( 228,i)*ze18 - pcm( 229,i)*ze19 - 
     &      pcm( 230,i)*ze20
c 
         ze22 = ze22 -
     &      pcm( 232,i)*ze1 - 
     &      pcm( 233,i)*ze2 - pcm( 234,i)*ze3 - pcm( 235,i)*ze4 - 
     &      pcm( 236,i)*ze5 - pcm( 237,i)*ze6 - pcm( 238,i)*ze7 - 
     &      pcm( 239,i)*ze8 - pcm( 240,i)*ze9 - pcm( 241,i)*ze10 - 
     &      pcm( 242,i)*ze11 - pcm( 243,i)*ze12 - pcm( 244,i)*ze13 - 
     &      pcm( 245,i)*ze14 - pcm( 246,i)*ze15 - pcm( 247,i)*ze16 - 
     &      pcm( 248,i)*ze17 - pcm( 249,i)*ze18 - pcm( 250,i)*ze19 - 
     &      pcm( 251,i)*ze20 - pcm( 252,i)*ze21
c 
         ze23 = ze23 -
     &      pcm( 254,i)*ze1 - 
     &      pcm( 255,i)*ze2 - pcm( 256,i)*ze3 - pcm( 257,i)*ze4 - 
     &      pcm( 258,i)*ze5 - pcm( 259,i)*ze6 - pcm( 260,i)*ze7 - 
     &      pcm( 261,i)*ze8 - pcm( 262,i)*ze9 - pcm( 263,i)*ze10 - 
     &      pcm( 264,i)*ze11 - pcm( 265,i)*ze12 - pcm( 266,i)*ze13 - 
     &      pcm( 267,i)*ze14 - pcm( 268,i)*ze15 - pcm( 269,i)*ze16 - 
     &      pcm( 270,i)*ze17 - pcm( 271,i)*ze18 - pcm( 272,i)*ze19 - 
     &      pcm( 273,i)*ze20 - pcm( 274,i)*ze21 - pcm( 275,i)*ze22
c 
         ze24 = ze24 -
     &      pcm( 277,i)*ze1 - 
     &      pcm( 278,i)*ze2 - pcm( 279,i)*ze3 - pcm( 280,i)*ze4 - 
     &      pcm( 281,i)*ze5 - pcm( 282,i)*ze6 - pcm( 283,i)*ze7 - 
     &      pcm( 284,i)*ze8 - pcm( 285,i)*ze9 - pcm( 286,i)*ze10 - 
     &      pcm( 287,i)*ze11 - pcm( 288,i)*ze12 - pcm( 289,i)*ze13 - 
     &      pcm( 290,i)*ze14 - pcm( 291,i)*ze15 - pcm( 292,i)*ze16 - 
     &      pcm( 293,i)*ze17 - pcm( 294,i)*ze18 - pcm( 295,i)*ze19 - 
     &      pcm( 296,i)*ze20 - pcm( 297,i)*ze21 - pcm( 298,i)*ze22 - 
     &      pcm( 299,i)*ze23
c
         z(edest(1,i)) = ze1
         z(edest(2,i)) = ze2
         z(edest(3,i)) = ze3
         z(edest(4,i)) = ze4
         z(edest(5,i)) = ze5
         z(edest(6,i)) = ze6
         z(edest(7,i)) = ze7
         z(edest(8,i)) = ze8
         z(edest(9,i)) = ze9
         z(edest(10,i)) = ze10
         z(edest(11,i)) = ze11
         z(edest(12,i)) = ze12
         z(edest(13,i)) = ze13
         z(edest(14,i)) = ze14
         z(edest(15,i)) = ze15
         z(edest(16,i)) = ze16
         z(edest(17,i)) = ze17
         z(edest(18,i)) = ze18
         z(edest(19,i)) = ze19
         z(edest(20,i)) = ze20
         z(edest(21,i)) = ze21
         z(edest(22,i)) = ze22
         z(edest(23,i)) = ze23
         z(edest(24,i)) = ze24
      end do
c$OMP END PARALLEL DO
c
      return
      end
