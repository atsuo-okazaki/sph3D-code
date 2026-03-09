      MODULE perform

      use idims
      implicit none
      save

      INTEGER(I4B) :: ninit1, ninit2, ninit3, ninit4, ninit5
      REAL(DP) :: tins, tins1, tins2, tginit1, tginit2, tginit, &
          tgload1, tgload2, tgload, tgcall11, tgcall12, tgcall1, tgsort1, &
          tgsort2, tgsort, tdens, tdens1, tdens2, tforce, tforce1, &
          tforce2, tggrav1, tggrav2, tggrav, &
          tgnei1, tgnei2, tgnei, tgcall21, tgcall22, tgcall2, &
          tgcall31, tgcall32, tgcall3, tgnmisc1, tgnmisc2, tgnmisc, &
          tgnstor1, tgnstor2, tgnstor, tgforpt1, tgforpt2, tgforpt, &
          ttest1, ttest2, ttest, tmtree1, tmtree2, tmtree, ttreef1, &
          ttreef2, ttreef, trevt1, trevt2, trevt

      END MODULE perform
