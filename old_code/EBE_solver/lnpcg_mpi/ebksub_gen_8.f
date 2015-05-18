c     ****************************************************************
c     *                                                              *
c     *                      subroutine ebksub_gen                   *
c     *                     **** version for MPI ****                *
c     *                       written by : bh                        *
c     *                                                              *
c     *             last modified : 11/27/09 rhd (no threads)        *
c     *                                                              *
c     *     this subroutine performs back substitution for a         *
c     *     block of non-conflicting, similar elements.              * 
c     *                                                              *
c     *     this version is the general case which handles any       *
c     *     element type.                                            *
c     *                                                              *
c     ****************************************************************
c     
      subroutine ebksub_gen ( span, totdof, edest, z, pcm, utsz )
      implicit integer (a-z)
c     
$add param_def
c     
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
c        cannot run thread parallel on MPI version. blocking is not
c        vectorized. threaded code here relies on no elements
c        in a block sharing a node for the correct result.
c
      do i = span, 1, -1
         do j = 1, totdof
            ze(j) = z(edest(j,i))
         end do
c
         do j = totdof - 1, 1, -1
            do k = totdof , j + 1, -1
               ze(j) = ze(j) - pcm( tops(k)+j-1,i) * ze(k)
            end do
         end do
c
         do j = 1, totdof
            z(edest(j,i)) = ze(j)
         end do
      end do
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ebksub8                      *
c     *                     **** version for MPI ****                *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/27/09 rhd (threads)     *
c     *                                                              *
c     *     this subroutine performs back substitution for a         *
c     *     block of non-conflicting, similar elements.              * 
c     *                                                              *
c     *     this version is coded specifically for 8 node elements.  *
c     *                                                              *
c     ****************************************************************
c     
c     
c     
      subroutine ebksub8 ( span, totdof, edest, z, pcm, utsz )
      implicit integer (a-z)
c     
$add param_def
c     
      dimension edest(totdof,*)
#sgl      real
#dbl      double precision
     &     pcm(utsz,*),z(*),
     &     ze1,ze2,ze3,ze4,ze5,ze6,ze7,ze8,ze9,
     &     ze10,ze11,ze12,ze13,ze14,ze15,ze16,ze17,ze18,ze19,
     &     ze20,ze21,ze22,ze23,ze24
c
c        cannot run thread parallel on MPI version. blocking is not
c        vectorized. threaded code here relies on no elements
c        in a block sharing a node for the correct result.
c
      do i = span, 1, -1
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
         ze23 = ze23 -
     &      pcm( 299,i)*ze24
c 
         ze22 = ze22 -
     &      pcm( 298,i)*ze24 - pcm( 275,i)*ze23
c 
         ze21 = ze21 -
     &      pcm( 297,i)*ze24 - pcm( 274,i)*ze23 - pcm( 252,i)*ze22
c 
         ze20 = ze20 -
     &      pcm( 296,i)*ze24 - pcm( 273,i)*ze23 - pcm( 251,i)*ze22 - 
     &      pcm( 230,i)*ze21
c 
         ze19 = ze19 -
     &      pcm( 295,i)*ze24 - pcm( 272,i)*ze23 - pcm( 250,i)*ze22 - 
     &      pcm( 229,i)*ze21 - pcm( 209,i)*ze20
c 
         ze18 = ze18 -
     &      pcm( 294,i)*ze24 - pcm( 271,i)*ze23 - pcm( 249,i)*ze22 - 
     &      pcm( 228,i)*ze21 - pcm( 208,i)*ze20 - pcm( 189,i)*ze19
c 
         ze17 = ze17 -
     &      pcm( 293,i)*ze24 - pcm( 270,i)*ze23 - pcm( 248,i)*ze22 - 
     &      pcm( 227,i)*ze21 - pcm( 207,i)*ze20 - pcm( 188,i)*ze19 - 
     &      pcm( 170,i)*ze18
c 
         ze16 = ze16 -
     &      pcm( 292,i)*ze24 - pcm( 269,i)*ze23 - pcm( 247,i)*ze22 - 
     &      pcm( 226,i)*ze21 - pcm( 206,i)*ze20 - pcm( 187,i)*ze19 - 
     &      pcm( 169,i)*ze18 - pcm( 152,i)*ze17
c 
         ze15 = ze15 -
     &      pcm( 291,i)*ze24 - pcm( 268,i)*ze23 - pcm( 246,i)*ze22 - 
     &      pcm( 225,i)*ze21 - pcm( 205,i)*ze20 - pcm( 186,i)*ze19 - 
     &      pcm( 168,i)*ze18 - pcm( 151,i)*ze17 - pcm( 135,i)*ze16
c 
         ze14 = ze14 -
     &      pcm( 290,i)*ze24 - pcm( 267,i)*ze23 - pcm( 245,i)*ze22 - 
     &      pcm( 224,i)*ze21 - pcm( 204,i)*ze20 - pcm( 185,i)*ze19 - 
     &      pcm( 167,i)*ze18 - pcm( 150,i)*ze17 - pcm( 134,i)*ze16 - 
     &      pcm( 119,i)*ze15
c 
         ze13 = ze13 -
     &      pcm( 289,i)*ze24 - pcm( 266,i)*ze23 - pcm( 244,i)*ze22 - 
     &      pcm( 223,i)*ze21 - pcm( 203,i)*ze20 - pcm( 184,i)*ze19 - 
     &      pcm( 166,i)*ze18 - pcm( 149,i)*ze17 - pcm( 133,i)*ze16 - 
     &      pcm( 118,i)*ze15 - pcm( 104,i)*ze14
c 
         ze12 = ze12 -
     &      pcm( 288,i)*ze24 - pcm( 265,i)*ze23 - pcm( 243,i)*ze22 - 
     &      pcm( 222,i)*ze21 - pcm( 202,i)*ze20 - pcm( 183,i)*ze19 - 
     &      pcm( 165,i)*ze18 - pcm( 148,i)*ze17 - pcm( 132,i)*ze16 - 
     &      pcm( 117,i)*ze15 - pcm( 103,i)*ze14 - pcm(  90,i)*ze13
c 
         ze11 = ze11 -
     &      pcm( 287,i)*ze24 - pcm( 264,i)*ze23 - pcm( 242,i)*ze22 - 
     &      pcm( 221,i)*ze21 - pcm( 201,i)*ze20 - pcm( 182,i)*ze19 - 
     &      pcm( 164,i)*ze18 - pcm( 147,i)*ze17 - pcm( 131,i)*ze16 - 
     &      pcm( 116,i)*ze15 - pcm( 102,i)*ze14 - pcm(  89,i)*ze13 - 
     &      pcm(  77,i)*ze12
c 
         ze10 = ze10 -
     &      pcm( 286,i)*ze24 - pcm( 263,i)*ze23 - pcm( 241,i)*ze22 - 
     &      pcm( 220,i)*ze21 - pcm( 200,i)*ze20 - pcm( 181,i)*ze19 - 
     &      pcm( 163,i)*ze18 - pcm( 146,i)*ze17 - pcm( 130,i)*ze16 - 
     &      pcm( 115,i)*ze15 - pcm( 101,i)*ze14 - pcm(  88,i)*ze13 - 
     &      pcm(  76,i)*ze12 - pcm(  65,i)*ze11
c 
         ze9 = ze9 -
     &      pcm( 285,i)*ze24 - pcm( 262,i)*ze23 - pcm( 240,i)*ze22 - 
     &      pcm( 219,i)*ze21 - pcm( 199,i)*ze20 - pcm( 180,i)*ze19 - 
     &      pcm( 162,i)*ze18 - pcm( 145,i)*ze17 - pcm( 129,i)*ze16 - 
     &      pcm( 114,i)*ze15 - pcm( 100,i)*ze14 - pcm(  87,i)*ze13 - 
     &      pcm(  75,i)*ze12 - pcm(  64,i)*ze11 - pcm(  54,i)*ze10
c 
         ze8 = ze8 -
     &      pcm( 284,i)*ze24 - pcm( 261,i)*ze23 - pcm( 239,i)*ze22 - 
     &      pcm( 218,i)*ze21 - pcm( 198,i)*ze20 - pcm( 179,i)*ze19 - 
     &      pcm( 161,i)*ze18 - pcm( 144,i)*ze17 - pcm( 128,i)*ze16 - 
     &      pcm( 113,i)*ze15 - pcm(  99,i)*ze14 - pcm(  86,i)*ze13 - 
     &      pcm(  74,i)*ze12 - pcm(  63,i)*ze11 - pcm(  53,i)*ze10 - 
     &      pcm(  44,i)*ze9
c 
         ze7 = ze7 -
     &      pcm( 283,i)*ze24 - pcm( 260,i)*ze23 - pcm( 238,i)*ze22 - 
     &      pcm( 217,i)*ze21 - pcm( 197,i)*ze20 - pcm( 178,i)*ze19 - 
     &      pcm( 160,i)*ze18 - pcm( 143,i)*ze17 - pcm( 127,i)*ze16 - 
     &      pcm( 112,i)*ze15 - pcm(  98,i)*ze14 - pcm(  85,i)*ze13 - 
     &      pcm(  73,i)*ze12 - pcm(  62,i)*ze11 - pcm(  52,i)*ze10 - 
     &      pcm(  43,i)*ze9 - pcm(  35,i)*ze8
c 
         ze6 = ze6 -
     &      pcm( 282,i)*ze24 - pcm( 259,i)*ze23 - pcm( 237,i)*ze22 - 
     &      pcm( 216,i)*ze21 - pcm( 196,i)*ze20 - pcm( 177,i)*ze19 - 
     &      pcm( 159,i)*ze18 - pcm( 142,i)*ze17 - pcm( 126,i)*ze16 - 
     &      pcm( 111,i)*ze15 - pcm(  97,i)*ze14 - pcm(  84,i)*ze13 - 
     &      pcm(  72,i)*ze12 - pcm(  61,i)*ze11 - pcm(  51,i)*ze10 - 
     &      pcm(  42,i)*ze9 - pcm(  34,i)*ze8 - pcm(  27,i)*ze7
c 
         ze5 = ze5 -
     &      pcm( 281,i)*ze24 - pcm( 258,i)*ze23 - pcm( 236,i)*ze22 - 
     &      pcm( 215,i)*ze21 - pcm( 195,i)*ze20 - pcm( 176,i)*ze19 - 
     &      pcm( 158,i)*ze18 - pcm( 141,i)*ze17 - pcm( 125,i)*ze16 - 
     &      pcm( 110,i)*ze15 - pcm(  96,i)*ze14 - pcm(  83,i)*ze13 - 
     &      pcm(  71,i)*ze12 - pcm(  60,i)*ze11 - pcm(  50,i)*ze10 - 
     &      pcm(  41,i)*ze9 - pcm(  33,i)*ze8 - pcm(  26,i)*ze7 - 
     &      pcm(  20,i)*ze6
c 
         ze4 = ze4 -
     &      pcm( 280,i)*ze24 - pcm( 257,i)*ze23 - pcm( 235,i)*ze22 - 
     &      pcm( 214,i)*ze21 - pcm( 194,i)*ze20 - pcm( 175,i)*ze19 - 
     &      pcm( 157,i)*ze18 - pcm( 140,i)*ze17 - pcm( 124,i)*ze16 - 
     &      pcm( 109,i)*ze15 - pcm(  95,i)*ze14 - pcm(  82,i)*ze13 - 
     &      pcm(  70,i)*ze12 - pcm(  59,i)*ze11 - pcm(  49,i)*ze10 - 
     &      pcm(  40,i)*ze9 - pcm(  32,i)*ze8 - pcm(  25,i)*ze7 - 
     &      pcm(  19,i)*ze6 - pcm(  14,i)*ze5
c 
         ze3 = ze3 -
     &      pcm( 279,i)*ze24 - pcm( 256,i)*ze23 - pcm( 234,i)*ze22 - 
     &      pcm( 213,i)*ze21 - pcm( 193,i)*ze20 - pcm( 174,i)*ze19 - 
     &      pcm( 156,i)*ze18 - pcm( 139,i)*ze17 - pcm( 123,i)*ze16 - 
     &      pcm( 108,i)*ze15 - pcm(  94,i)*ze14 - pcm(  81,i)*ze13 - 
     &      pcm(  69,i)*ze12 - pcm(  58,i)*ze11 - pcm(  48,i)*ze10 - 
     &      pcm(  39,i)*ze9 - pcm(  31,i)*ze8 - pcm(  24,i)*ze7 - 
     &      pcm(  18,i)*ze6 - pcm(  13,i)*ze5 - pcm(   9,i)*ze4
c 
         ze2 = ze2 -
     &      pcm( 278,i)*ze24 - pcm( 255,i)*ze23 - pcm( 233,i)*ze22 - 
     &      pcm( 212,i)*ze21 - pcm( 192,i)*ze20 - pcm( 173,i)*ze19 - 
     &      pcm( 155,i)*ze18 - pcm( 138,i)*ze17 - pcm( 122,i)*ze16 - 
     &      pcm( 107,i)*ze15 - pcm(  93,i)*ze14 - pcm(  80,i)*ze13 - 
     &      pcm(  68,i)*ze12 - pcm(  57,i)*ze11 - pcm(  47,i)*ze10 - 
     &      pcm(  38,i)*ze9 - pcm(  30,i)*ze8 - pcm(  23,i)*ze7 - 
     &      pcm(  17,i)*ze6 - pcm(  12,i)*ze5 - pcm(   8,i)*ze4 - 
     &      pcm(   5,i)*ze3
c 
         ze1 = ze1 -
     &      pcm( 277,i)*ze24 - pcm( 254,i)*ze23 - pcm( 232,i)*ze22 - 
     &      pcm( 211,i)*ze21 - pcm( 191,i)*ze20 - pcm( 172,i)*ze19 - 
     &      pcm( 154,i)*ze18 - pcm( 137,i)*ze17 - pcm( 121,i)*ze16 - 
     &      pcm( 106,i)*ze15 - pcm(  92,i)*ze14 - pcm(  79,i)*ze13 - 
     &      pcm(  67,i)*ze12 - pcm(  56,i)*ze11 - pcm(  46,i)*ze10 - 
     &      pcm(  37,i)*ze9 - pcm(  29,i)*ze8 - pcm(  22,i)*ze7 - 
     &      pcm(  16,i)*ze6 - pcm(  11,i)*ze5 - pcm(   7,i)*ze4 - 
     &      pcm(   4,i)*ze3 - pcm(   2,i)*ze2
         
c
         z(edest(1,i))  = ze1 
         z(edest(2,i))  = ze2 
         z(edest(3,i))  = ze3 
         z(edest(4,i))  = ze4 
         z(edest(5,i))  = ze5 
         z(edest(6,i))  = ze6 
         z(edest(7,i))  = ze7 
         z(edest(8,i))  = ze8 
         z(edest(9,i))  = ze9 
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
c
      return
      end
