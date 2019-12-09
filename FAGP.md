FAGP
================
Analine Aguayo
12/4/2019

``` r
library(bio3d)
read.fasta("Sequences.fasta.txt", rm.dup=FALSE,)
```

    ##                                1        .         .         .         .         50 
    ## [Truncated_Name:1]Human        MASQPNSSAKKKEEKGKNIQVVVRCRPFNLAERKASAHSIVECDPVRKEV
    ## [Truncated_Name:2]Three-spin   FXXXXXXXFGPESRDQDEVYQILERGSAKRRTASTLMNAYSSRSHSVFSV
    ## [Truncated_Name:3]Clownfish    MASHNQSGAKREEKGRNIQVVVRCRPFNTVERKSSYGVIDCDQSRKEVMV
    ## [Truncated_Name:4]Gilt-head_   MASHNLSGVKREEKGRNIQVVVRCRPFNTMERKSSYGVIDCDSNRKEVMV
    ## [Truncated_Name:5]House_cat    MASQPNSSAKKKEEKGKNIQVVVRCRPFNLAERKASAHSVVECDHVRKEV
    ## [Truncated_Name:6]Mouse        MASQPSSLKKKEEKGRNIQVVVRCRPFNLAERKANAHSVVECDHARKEVS
    ## [Truncated_Name:7]Yellow_per   MSSHNLSGLKRDEKGRNIQVVVRCRPFNTMERKSSYGVMDCEPNRKEVSV
    ## [Truncated_Name:8]Zander_fis   MASHNLSGLKRDEKGRNIQVVVRCRPFNTMERKSSYGVMDCEPSRKEVSV
    ## [Truncated_Name:9]Zebrafish    MASSQVPAAKKDEKGRNIQVVVRCRPFNTVERKSGSHTVVECDQNRKEVI
    ##                                                     ^                             
    ##                                1        .         .         .         .         50 
    ## 
    ##                               51        .         .         .         .         100 
    ## [Truncated_Name:1]Human        SVRTGGLADKSSRKTYTFDMVFGASTKQIDVYRSVVCPILDEVIMGYNCT
    ## [Truncated_Name:2]Three-spin   TIHMKEITMDGEELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQS
    ## [Truncated_Name:3]Clownfish    KTGGMNDKASRKTYTFDMVFGPAAKQIDVYRSVVCPILDEVIMGYNCTVF
    ## [Truncated_Name:4]Gilt-head_   KTGGMNDKASRKTYTFDMVFGPAAKQIEVYRNVVCPILDEVIMGYNCTVF
    ## [Truncated_Name:5]House_cat    SVRTGGLADKSSRKTYTFDMVFGASTKQIDVYRSVVCPILDEVIMGYNCT
    ## [Truncated_Name:6]Mouse        VRTAGLTDKTSKKTYTFDMVFGASTKQIDVYRSVVCPILDEVIMGYNCTI
    ## [Truncated_Name:7]Yellow_per   RTGGMNDKAARKTYTFDMVFGQAAKQIDVYRSVVCPILDEVIMGYNCTVF
    ## [Truncated_Name:8]Zander_fis   RTGGMNDKAARKTYTFDMVFGQAAKQIDVYRSVVCPILDEVIMGYNCTVF
    ## [Truncated_Name:9]Zebrafish    MRTGGATDKAARKTYTFDMVFGPSAKQIEVYRSVVCPILDEVIMGYNCTI
    ##                                                                                   
    ##                               51        .         .         .         .         100 
    ## 
    ##                              101        .         .         .         .         150 
    ## [Truncated_Name:1]Human        IFAYGQTGTGKTFTMEGERSPNEEYTWEEDPLAGIIPRTLHQIFEKLTDN
    ## [Truncated_Name:2]Three-spin   LLTLGRVITALVEKRPHVPYRESKLTRILQDSLGGRTKTSIIATVSPSSS
    ## [Truncated_Name:3]Clownfish    AYGQTGTGKTFTMEGERSPDGEFTWEEDPLAGIIPRTLHQIFEKLSENGT
    ## [Truncated_Name:4]Gilt-head_   AYGQTGTGKTFTMEGERSPNEQFTWEEDPLAGIIPRTLHQIFEKLSENGT
    ## [Truncated_Name:5]House_cat    IFAYGQTGTGKTFTMEGERSPNEEYTWEEDPLAGIIPRTLHQIFEKLTDN
    ## [Truncated_Name:6]Mouse        FAYGQTGTGKTFTMEGERSPNEVYTWEEDPLAGIIPRTLHQIFEKLTDNG
    ## [Truncated_Name:7]Yellow_per   AYGQTGTGKTFTMEGERSPDGEFTWEEDPLAGIIPRTLHQIFEKLSENGT
    ## [Truncated_Name:8]Zander_fis   AYGQTGTGKTFTMEGERSPDGEFTWEEDPLAGIIPRTLHQIFEKLSENGT
    ## [Truncated_Name:9]Zebrafish    FAYGQTGTGKTFTMEGDRSPNEEFTWEEDPLAGIIPRTLHQIFEKLSNNG
    ##                                                                                   
    ##                              101        .         .         .         .         150 
    ## 
    ##                              151        .         .         .         .         200 
    ## [Truncated_Name:1]Human        GTEFSVKVSLLEIYNEELFDLLNPSSDVSERLQMFDDPRNKRGVIIKGLE
    ## [Truncated_Name:2]Three-spin   NLEETLSTLEYASRAKNIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAA
    ## [Truncated_Name:3]Clownfish    EFSVKVSLLEIYNEELFDLLSPTEDVNERLQLFDDPRNKRGVVVKGLEEV
    ## [Truncated_Name:4]Gilt-head_   EFSVKVSLLEIYNEELFDLLSPTEDVNERLQLFDDPRNKRGVVVKGLEEV
    ## [Truncated_Name:5]House_cat    GTEFSVKVSLLEIYNEELFDLLNPSSDVSERLQMFDDPRNKRGVIIKGLE
    ## [Truncated_Name:6]Mouse        TEFSVKVSLLEIYNEELFDLLSPSSDVSERLQMFDDPRNKRGVIIKGLEE
    ## [Truncated_Name:7]Yellow_per   EFSVKVSLLEIYNEELFDLLSPSDDVNERLQLFDDPRNKRGVVVKGLEEV
    ## [Truncated_Name:8]Zander_fis   EFSVKVSLLEIYNEELFDLLSPSDDVNERLQLFDDPRNKRGVVVKGLEEV
    ## [Truncated_Name:9]Zebrafish    TEFSVKVSLLEIYNEELFDLLSPAPDVTERLQLFDDPRNKRGVTIKGLEE
    ##                                                                                   
    ##                              151        .         .         .         .         200 
    ## 
    ##                              201        .         .         .         .         250 
    ## [Truncated_Name:1]Human        EITVHNKDEVYQILEKGAAKRTTAATLMNAYSSRSHSVFSVTIHMKETTI
    ## [Truncated_Name:2]Three-spin   TRDKNGIYLSAENYESMMGQITSHEVHTVEYSDRIAAMEEEIKKVTELFV
    ## [Truncated_Name:3]Clownfish    TVHNKDEVYQILERGAAKRRTASTLMNAYSSRSHSVFSVTIHMKEITVDG
    ## [Truncated_Name:4]Gilt-head_   TVHNKDEVYQILERGSAKRRTASTLMNAYSSRSHSVFSVTIHMKEITLDG
    ## [Truncated_Name:5]House_cat    EITVHNKDEVYQILEKGAAKRTTAATLMNAYSSRSHSVFSVTIHMKETTI
    ## [Truncated_Name:6]Mouse        ITVHNKDEVYQILEKGAAKRTTAATLMNAYSSRSHSVFSVTIHMKETTID
    ## [Truncated_Name:7]Yellow_per   TVHNKDEVYQILERGSAKRRTASTLMNAYSSRSHSVFSVTIHMKEITMDG
    ## [Truncated_Name:8]Zander_fis   TVHNKDEVYQILERGSAKRRTASTLMNAYSSRSHSVFSVTIHMKEITMDG
    ## [Truncated_Name:9]Zebrafish    ITVHNKNEVYQILERGAAKRKTASTLMNAYSSRSHSVFSVTIHMKEITLD
    ##                                                                                   
    ##                              201        .         .         .         .         250 
    ## 
    ##                              251        .         .         .         .         300 
    ## [Truncated_Name:1]Human        DGEELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVIT
    ## [Truncated_Name:2]Three-spin   DSKTRLELCAVDLDEKQQRLEETSRDLQHTKEKLMEXEFVCSELTLVQES
    ## [Truncated_Name:3]Clownfish    EELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVITAL
    ## [Truncated_Name:4]Gilt-head_   EELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVITAL
    ## [Truncated_Name:5]House_cat    DGEELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVIT
    ## [Truncated_Name:6]Mouse        GEELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVITA
    ## [Truncated_Name:7]Yellow_per   EELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVITAL
    ## [Truncated_Name:8]Zander_fis   EELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVITAL
    ## [Truncated_Name:9]Zebrafish    GEELVKIGKLNLVDLAGSENIGRSGAVDKRAREAGNINQSLLTLGRVIKA
    ##                                                                                   
    ##                              251        .         .         .         .         300 
    ## 
    ##                              301        .         .         .         .         350 
    ## [Truncated_Name:1]Human        ALVERTPHVPYRESKLTRILQDSLGGRTRTSIIATISPASLNLEETLSTL
    ## [Truncated_Name:2]Three-spin   LYDTAGRLLSTVDASTGDVCGLPGQLDRXKXVEQHYSGVQQSSLSAWX--
    ## [Truncated_Name:3]Clownfish    VEKRPHIPYRESKLTRILQDSLGGRTKTSIIATVSPSSSNLEETLSTLEY
    ## [Truncated_Name:4]Gilt-head_   VEKRPHVPYRESKLTRILQDSLGGRTKTSIIATVSPSSSNLEETLSTLEY
    ## [Truncated_Name:5]House_cat    ALVERTPHVPYRESKLTRILQDSLGGRTRTSIIATISPASLNLEETLSTL
    ## [Truncated_Name:6]Mouse        LVERTPHIPYRESKLTRILQDSLGGRTRTSIIATISPASFNLEETLSTLE
    ## [Truncated_Name:7]Yellow_per   VEKRPHIPYRESKLTRILQDSLGGRTKTSIIATVSPSSSNLEETLSTLEY
    ## [Truncated_Name:8]Zander_fis   VEKRPHIPYRESKLTRILQDSLGGRTKTSIIATVSPSSSNLEETLSTLEY
    ## [Truncated_Name:9]Zebrafish    LVERGPHVPYRESKLTRILQDSLGGRTKTSIIATVSPASINLEETLSTLD
    ##                                                                                   
    ##                              301        .         .         .         .         350 
    ## 
    ##                              351        .         .         .         .         400 
    ## [Truncated_Name:1]Human        EYAHRAKNILNKPEVNQKLTKKALIKEYTEEIERLKRDLAAAREKNGVYI
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    ASRAKNIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAATRDKNGVYLSA
    ## [Truncated_Name:4]Gilt-head_   ASRAKNIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAATRDKNGVYLSA
    ## [Truncated_Name:5]House_cat    EYAHRAKNILNKPEVNQKLTKRALIKEYTEEIERLKRDLAAAREKNGVYI
    ## [Truncated_Name:6]Mouse        YAHRAKNIMNKPEVNQKLTKKALIKEYTEEIERLKRDLAAAREKNGVYIS
    ## [Truncated_Name:7]Yellow_per   ASRAKNIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAATRDKNGVYLSA
    ## [Truncated_Name:8]Zander_fis   ASRAKNIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAATRDKNGVYLSA
    ## [Truncated_Name:9]Zebrafish    YANRAKSIMNKPEVNQKLTKRTLIKEYTEEIERLKRDLAATRDKHGVYLS
    ##                                                                                   
    ##                              351        .         .         .         .         400 
    ## 
    ##                              401        .         .         .         .         450 
    ## [Truncated_Name:1]Human        SEENFRVMSGKLTVQEEQIVELIEKIGAVEEELNRVTELFMDNKNELDQC
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    ENYESMVGQITSHEEQIVEYTDKIAAMEEEIKKVTELFEESKTRLEQCTV
    ## [Truncated_Name:4]Gilt-head_   ENYESMMGQITSQEEQTTECTERIAAMEEELRKVTELFEDSKSKLEQCAG
    ## [Truncated_Name:5]House_cat    SEENFRAMSGKLTVQEEQIVELIEKIGAVEEELSRVTELFMDNKNELDQC
    ## [Truncated_Name:6]Mouse        EESFRAMNGKVTVQEEQIVELVEKIAVLEEELSKATELFMDSKNELDQCK
    ## [Truncated_Name:7]Yellow_per   ENYESMMGQITSHEEHTSEYTDRIAAMEEEIKKVTELFMDSKTRLEQCTV
    ## [Truncated_Name:8]Zander_fis   ENYESMMGQITSHEEHTSEYTDRIAAMEEEIKKVTELFTDSKTRLEQCTV
    ## [Truncated_Name:9]Zebrafish    VDNYETLNGKIVSQEEQITEYTERIAAMEEELKKIIDLFTDSKQKLEQCT
    ##                                                                                   
    ##                              401        .         .         .         .         450 
    ## 
    ##                              451        .         .         .         .         500 
    ## [Truncated_Name:1]Human        KSDLQNKTQELETTQKHLQETKLQLVKEEYITSALESTEEKLHDAASKLL
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    DLDDKQQRLEETSRDLQQTKEKLSEEEFICSELTSAQETLYNAAGQLLST
    ## [Truncated_Name:4]Gilt-head_   ELDEKQQRLEETSQDLQQTKEKLSEEEFISSELSSVQTTLYSTAGQLLST
    ## [Truncated_Name:5]House_cat    KSDLQNKTQELETTQKHLQETKLQLVEEEYITSALESTEEKLHDAASRLL
    ## [Truncated_Name:6]Mouse        SDLQTKTQELETTQKHLQETKLQLVKEEYVSSALERTEKTLHDTASKLLN
    ## [Truncated_Name:7]Yellow_per   DLDEKQQILEETSKDLQQTKEKLSQEEFVCSELAVVQETLYNTAGQLLST
    ## [Truncated_Name:8]Zander_fis   DLDQKQQMLEETSKDLQQTKEKLSQEEFVCSELTVVQETLYNTAGQLLST
    ## [Truncated_Name:9]Zebrafish    EDLQDKNQRLEEAHKDLSETRHRLNQEEFISTQLQTNESHLYNTADQLLS
    ##                                                                                   
    ##                              451        .         .         .         .         500 
    ## 
    ##                              501        .         .         .         .         550 
    ## [Truncated_Name:1]Human        NTVEETTKDVSGLHSKLDRKKAVDQHNAEAQDIFGKNLNSLFNNMEELIK
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    VDASTSDVSGLHDKLDRKKKVEQHNSQMQQTFAQRMDRAFSSLQCSIQQH
    ## [Truncated_Name:4]Gilt-head_   ADASTSDVTGLHDKLDRKKKVEQHNSQIQQSFAQRMDGAFSDMQRRVQQQ
    ## [Truncated_Name:5]House_cat    TTVEETTKDVSGLHSKLDRKKAIDQHNAEAQDVFGKNLNGLFNRMEELIK
    ## [Truncated_Name:6]Mouse        TVKETTRAVSGLHSKLDRKRAIDEHNAEAQESFGKNLNSLFNNMEELIKD
    ## [Truncated_Name:7]Yellow_per   VDASTSDVMGLHDKLDRKKKVEQHNSQIQQSFSQRMEGALSNMQRCVQQH
    ## [Truncated_Name:8]Zander_fis   VDASTSDVMGLHDKLDRKKKVEQHNSQIQQSFSQRMEGALSNMQRCVQQH
    ## [Truncated_Name:9]Zebrafish    TAEASTQDVGGLHAKLQRKKDVELHNSKVQESFSQCMENCYNSMQTSLKE
    ##                                                                                   
    ##                              501        .         .         .         .         550 
    ## 
    ##                              551        .         .         .         .         600 
    ## [Truncated_Name:1]Human        DGSSKQKAMLEVHKTLFGNLLSSSVSALDTITTVALGSLTSIPENVSTHV
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    GDKHRDMLSGYSQAIDGLLVKNEAALKGALTTVESFVDGVGPLVAEGVSR
    ## [Truncated_Name:4]Gilt-head_   GTKHHDMLSSYSQAVDGLLVMNEAVLKGALTSVESLVGGVGQLVAEGVAR
    ## [Truncated_Name:5]House_cat    DSSSKQKVMLEVHKTLFGNLLSSSVSALDTITTTALGSLTSVPENVSTHV
    ## [Truncated_Name:6]Mouse        GSAKQKAMLDVHKTLFGNLMSSSVSALDTITTTALESLVSIPENVSARVS
    ## [Truncated_Name:7]Yellow_per   GDKHHDMLSSYSQAIDGVRVMNEAALKGALTTVEAFVGGIGTLVADGVSR
    ## [Truncated_Name:8]Zander_fis   GDKHHDMLSSYSQAIDGVRVMNEAALKGALTTVESFVGGVGTLVADGVSR
    ## [Truncated_Name:9]Zebrafish    QSQKHAAMIDYYRSSVGELLNTNGKVFKETLSAVCESYSSIKGAVGEGVE
    ##                                                                                   
    ##                              551        .         .         .         .         600 
    ## 
    ##                              601        .         .         .         .         650 
    ## [Truncated_Name:1]Human        SQIFNMILKEQSLAAESKTVLQELINVLKTDLLSSLEMILSPTVVSILKI
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    CQERVQQQEALCLQDKEALLQLLEEHRQEVEELLVARTLMGLSAVKELSD
    ## [Truncated_Name:4]Gilt-head_   CRDKVQQQEALCLQDKESMLQLLEEHRQDMEEVLVTRTLMGLSAAKELND
    ## [Truncated_Name:5]House_cat    SQISNMILKEQSLAAESKTVLQTLIDVLKTDLLSSLETTLSPTVVSILKI
    ## [Truncated_Name:6]Mouse        QISDMILEEQSLAAQSKSVLQGLIDELVTDLFTSLKTIVAPSVVSILNIN
    ## [Truncated_Name:7]Yellow_per   CRDRVQQQEALCLQDKESLLQLLEEHRQDMEDVLVARTLMGLSAVKELND
    ## [Truncated_Name:8]Zander_fis   CRDRVQQQEALCLQDKESLLQLLEEHRQDMEDVLVARTLMGLSAVKELND
    ## [Truncated_Name:9]Zebrafish    RCKEQVLNQEKLSQDAQNSILEILDEHKQHLEEVLVAQAVPGIRSVMSMN
    ##                                                                                   
    ##                              601        .         .         .         .         650 
    ## 
    ##                              651        .         .         .         .         700 
    ## [Truncated_Name:1]Human        NSQLKHIFKTSLTVADKIEDQKKELDGFLSILCNNLHELQENTICSLVES
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    TLRASVETQRVLADKVEAMKEMGVFFSGLVRELAGLRQASVQGLSDLQAE
    ## [Truncated_Name:4]Gilt-head_   SLRAAVDTQRALADKVEAMKEVGVFFSGLVQELAGLREAAVKGLSSLQAE
    ## [Truncated_Name:5]House_cat    NSQLKHIFKTSLTVASKIEDQKKEMDGFLSMLCNSLHELRENTISSLAES
    ## [Truncated_Name:6]Mouse        KQLQHIFRASSTVAEKVEDQKREIDSFLSILCNNLHELRENTVSSLVESQ
    ## [Truncated_Name:7]Yellow_per   TLRATVETQRALADKVEAMKEVGVVFGGLVRELAGLREEAVRGLSSLQAE
    ## [Truncated_Name:8]Zander_fis   TLRATVETQRALADKVEAMKELGVSFGGLVRELAGLREEAVQGLSSLQAE
    ## [Truncated_Name:9]Zebrafish    DNLKQTLHKYHNLAEQMQGVKADMMTFFDAYTESLASMRECALQGFNTLR
    ##                                                                                   
    ##                              651        .         .         .         .         700 
    ## 
    ##                              701        .         .         .         .         750 
    ## [Truncated_Name:1]Human        QKQCGNLTEDLKTIKQTHSQELCKLMNLWTERFCALEEKCENIQKPLSSV
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    HDKLEDEIRQAQERHQTGMKQTIQCLQDQLNLLNMEAQKDYMDLRSASKA
    ## [Truncated_Name:4]Gilt-head_   QDELEDKIRQAQERHHTGMKQTIQCLQDQLNLLSMEAEKDFTDLRSASRS
    ## [Truncated_Name:5]House_cat    QKLCENLTEDLKTIKKIHSQELDQLISLWAERFCALEEKCENIQKPLSSV
    ## [Truncated_Name:6]Mouse        KLCGDLTEDLKTIKETHSQELCQLSSSWAERFCALEKKYENIQKPLNSIQ
    ## [Truncated_Name:7]Yellow_per   HGKLEDEIRRAQDRHQTGMKQTMQLLQEQLNLLTMESQKDYSDLRSASKA
    ## [Truncated_Name:8]Zander_fis   HGKLEDEIRQAQDRHQTGMKQTMQLLQEQLNLLTMEAQKDYADLRSASTA
    ## [Truncated_Name:9]Zebrafish    AEHDKLKQQISQAGNSHQVRVAELVQCLQNQMNLLAVDTQNDFEGLSQAA
    ##                                                                                   
    ##                              701        .         .         .         .         750 
    ## 
    ##                              751        .         .         .         .         800 
    ## [Truncated_Name:1]Human        QENIQQKSKDIVNKMTFHSQKFCADSDGFSQELRNFNQEGTKLVEESVKH
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    VQKPLQSLQENISSGCSSVEGQASAQADLLSSTSSSIASSLRLTADESRQ
    ## [Truncated_Name:4]Gilt-head_   LQKPLQVLQENISSGCSSVERQASVQADLLSSTSSSLASSLRQTADESRQ
    ## [Truncated_Name:5]House_cat    QENTEQKSKDIISKTTSYNTQFCADCDGLSQELRHFNQEGTKLVEESVKH
    ## [Truncated_Name:6]Mouse        ENTELRSTDIINKTTVHSKKILAESDGLLQELRHFNQEGTQLVEESVGHC
    ## [Truncated_Name:7]Yellow_per   LQKPLQTVQENISSGCSAVEHQASAQADLLSSTSAFLASSLRLSADESQQ
    ## [Truncated_Name:8]Zander_fis   LQKPLQTVQENISSGCSTVERQASAQADLLSSTSASLASSLRLSADESQQ
    ## [Truncated_Name:9]Zebrafish    SAQIPPLETLQSSIESKCTVAEEQAVSVRARLGSSVHGVISEMNSVTKEG
    ##                                                                                   
    ##                              751        .         .         .         .         800 
    ## 
    ##                              801        .         .         .         .         850 
    ## [Truncated_Name:1]Human        SDKLNGNLEKISQETEQRCESLNTRTVYFSEQWVSSLNEREQELHNLLEV
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    MLEEMSGCCSHLHNSVSELVQRDLQWSSRAREQAETQAQEHLTLMGKVST
    ## [Truncated_Name:4]Gilt-head_   TLEEMTGCCAHLHSSVSGLVERDLQWSSRVREHAETQAQENLSLMESIST
    ## [Truncated_Name:5]House_cat    CDKLNSSLELVSQETEQRCEALNKSTLSFSEQWVSCLNIREEELQNLLEV
    ## [Truncated_Name:6]Mouse        SSLNSNLETVSQEITQKCGTLNTSTVHFSDQWASCLSKRKEELENLMEFV
    ## [Truncated_Name:7]Yellow_per   ALDEMTGCCSHLHSSVSGLVQRDLQWTSRAREHAETQAQEHLSLSGKICT
    ## [Truncated_Name:8]Zander_fis   ALDEMTGCCSHLHGSVSGLVQRDLQWTSRAREHAETQAQEHLSLSGKIST
    ## [Truncated_Name:9]Zebrafish    ERALEECAGYCGHLQTSLDSLAESGLKWCDEAKGLTESKAQEQLKLIRQT
    ##                                                                                   
    ##                              801        .         .         .         .         850 
    ## 
    ##                              851        .         .         .         .         900 
    ## [Truncated_Name:1]Human        VSQCCEASSSDITEKSDGRKAAHEKQHNIFLDQMTIDEDKLIAQNLELNE
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    EAQTLRQSVELRCTEQLRAAEEDISRQQEEVKQALMKVRNQTSLDRTVLD
    ## [Truncated_Name:4]Gilt-head_   EVQTLHQNVQTRCTEQLRAAGQELSGRQEEVKRSLTAVQDQTAENRTVLD
    ## [Truncated_Name:5]House_cat    VNQGCEASSSEITEKLNEHKAANENQRNMFFDQITTDEEKLIAQSLEFNE
    ## [Truncated_Name:6]Mouse        NGCCKASSSEITKKVREQSAAVANQHSSFVAQMTSDEESCKAGSLELDKT
    ## [Truncated_Name:7]Yellow_per   EAQTLCQSVETRCAEQLRAAEEELSSRQAEVERSLVAVQDQTSADRTVLD
    ## [Truncated_Name:8]Zander_fis   EAQTLCQSVATRCTEQLRAAEEELSSRQAEVTRSLVAVQNQTSADRTVLN
    ## [Truncated_Name:9]Zebrafish    DTAVQDLLKSVEEKGEKAVQDCEARLGQMQQEMEATLGRVEMQTSKDEAT
    ##                                                                                   
    ##                              851        .         .         .         .         900 
    ## 
    ##                              901        .         .         .         .         950 
    ## [Truncated_Name:1]Human        TIKIGLTKLNCFLEQDLKLDIPTGTTPQRKSYLYPSTLVRTEPREHLLDQ
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    EQQAELQSSVETSQQLVHGFLQDELQQDVPTGATPQRREFVYPQQLVKSR
    ## [Truncated_Name:4]Gilt-head_   QQQAELQDHMETCQQLVHGFLQEELQQDVPTGATPQRREFVYPRKLMKLQ
    ## [Truncated_Name:5]House_cat    TIKIGLTKLNCFLQQDLKLDIPTGTTPQRKNYLYPSTLVRTEPREQLLDQ
    ## [Truncated_Name:6]Mouse        IKTGLTKLNCFLKQDLKLDIPTGMTPERKKYLYPTTLVRTEPREQLLDQL
    ## [Truncated_Name:7]Yellow_per   QQQAELRDHVETSQELVQGFLQEELQQDLPTGATPQRREFVYPRQLEKSR
    ## [Truncated_Name:8]Zander_fis   QQQAELRDHVETSQDLVQGFLQEELQQDLPTGATPQRREFVYPRQLEKSR
    ## [Truncated_Name:9]Zebrafish    LQEHRETLSSINTQALDTVHNFISSELRQDLPTGTTPQRKEYMYPRVLSR
    ##                                                                                   
    ##                              901        .         .         .         .         950 
    ## 
    ##                              951        .         .         .         .         1000 
    ## [Truncated_Name:1]Human        LKRKQPELLMMLNCSENNKEETIPDVDVEEAVLGQYTEEPLSQEPSVDAG
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    SRSELLESLRRQQEELQAAMEEEEEEEEEELEEDDKHSQDSVEDDMSTCN
    ## [Truncated_Name:4]Gilt-head_   SRSELLESLRRQQEELQAAMEEEELDEEEEEDKHSQVDHDSLEDELSTCN
    ## [Truncated_Name:5]House_cat    LKRKQPELLMMLNCSENNKEEISQDLDEEKSVLGHSIEEPPSQEPSIDAS
    ## [Truncated_Name:6]Mouse        QKKQPPMMLNSSEASKETSQDMDEEREALEQCTEELVSPETTEHPSADCS
    ## [Truncated_Name:7]Yellow_per   SRSELLESLRRRQEELQAAMEKEHHQEEEEEEVDHDSLEDEISTCNESVA
    ## [Truncated_Name:8]Zander_fis   SRSELLESLRRRQEELQVAMEKEHHQEEEEEVDHDSLEDELSTCNESLAT
    ## [Truncated_Name:9]Zebrafish    PRSREELEEEFRAQQEQLQSELKPCEIVMEVEEEKPVDQDSLEDDVSVSS
    ##                                                                                   
    ##                              951        .         .         .         .         1000 
    ## 
    ##                             1001        .         .         .         .         1050 
    ## [Truncated_Name:1]Human        VDCSSIGGVPFFQHKKSHGKDKENRGINTLERSKVEETTEHLVTKSRLPL
    ## [Truncated_Name:2]Three-spin   --------------------------------------------------
    ## [Truncated_Name:3]Clownfish    ESLATEPSFIDENLVFNESKRVPFFKKKKGGKKESKTSSRCKASENTSTP
    ## [Truncated_Name:4]Gilt-head_   ESVATEPSFIDENLVFNESKRVPFFKKKKGGKKDTKIPARPKASDNDAAS
    ## [Truncated_Name:5]House_cat    VDCSSSGGIPFFQHKKITWKRQRKQGH-----------------------
    ## [Truncated_Name:6]Mouse        SSRGLPFFQRKKPHGKDKENRGLNPVEKYKVEEASDLSISKSRLPLHTSI
    ## [Truncated_Name:7]Yellow_per   TEPSFIDENLVFNESKRVPFFKQKKSSKKEAKLLAKYKASENEASTPQKS
    ## [Truncated_Name:8]Zander_fis   EPSFIDENLVFNESKRVPFFKQKKSSKKEAKMLAKYKASENDATTPQKSK
    ## [Truncated_Name:9]Zebrafish    DGNNTEQSCSDENLICYENGRIPFFKKKSKKENGSKSLNRSKVENDSMST
    ##                                                                                   
    ##                             1001        .         .         .         .         1050 
    ## 
    ##                             1051        .  1063 
    ## [Truncated_Name:1]Human        RAQINL-------
    ## [Truncated_Name:2]Three-spin   -------------
    ## [Truncated_Name:3]Clownfish    QKSRLPLRCQN--
    ## [Truncated_Name:4]Gilt-head_   TPQKSRLPLRCQN
    ## [Truncated_Name:5]House_cat    -------------
    ## [Truncated_Name:6]Mouse        NL-----------
    ## [Truncated_Name:7]Yellow_per   KLPLRCQN-----
    ## [Truncated_Name:8]Zander_fis   LPLRCQN------
    ## [Truncated_Name:9]Zebrafish    PPRSKLPLRCQN-
    ##                                              
    ##                             1051        .  1063 
    ## 
    ## Call:
    ##   read.fasta(file = "Sequences.fasta.txt", rm.dup = FALSE)
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   9 sequence rows; 1063 position columns (348 non-gap, 715 gap) 
    ## 
    ## + attr: id, ali, call

``` r
x <- read.fasta("Sequences.fasta.txt", rm.dup=FALSE)
```

``` r
head(x$ali[,1:348])
```

    ##                               [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
    ## Human                         "M"  "A"  "S"  "Q"  "P"  "N"  "S"  "S"  "A" 
    ## Three-spined_stickleback_fish "F"  "X"  "X"  "X"  "X"  "X"  "X"  "X"  "F" 
    ## Clownfish                     "M"  "A"  "S"  "H"  "N"  "Q"  "S"  "G"  "A" 
    ## Gilt-head_bream_fish          "M"  "A"  "S"  "H"  "N"  "L"  "S"  "G"  "V" 
    ## House_cat                     "M"  "A"  "S"  "Q"  "P"  "N"  "S"  "S"  "A" 
    ## Mouse                         "M"  "A"  "S"  "Q"  "P"  "S"  "S"  "L"  "K" 
    ##                               [,10] [,11] [,12] [,13] [,14] [,15] [,16]
    ## Human                         "K"   "K"   "K"   "E"   "E"   "K"   "G"  
    ## Three-spined_stickleback_fish "G"   "P"   "E"   "S"   "R"   "D"   "Q"  
    ## Clownfish                     "K"   "R"   "E"   "E"   "K"   "G"   "R"  
    ## Gilt-head_bream_fish          "K"   "R"   "E"   "E"   "K"   "G"   "R"  
    ## House_cat                     "K"   "K"   "K"   "E"   "E"   "K"   "G"  
    ## Mouse                         "K"   "K"   "E"   "E"   "K"   "G"   "R"  
    ##                               [,17] [,18] [,19] [,20] [,21] [,22] [,23]
    ## Human                         "K"   "N"   "I"   "Q"   "V"   "V"   "V"  
    ## Three-spined_stickleback_fish "D"   "E"   "V"   "Y"   "Q"   "I"   "L"  
    ## Clownfish                     "N"   "I"   "Q"   "V"   "V"   "V"   "R"  
    ## Gilt-head_bream_fish          "N"   "I"   "Q"   "V"   "V"   "V"   "R"  
    ## House_cat                     "K"   "N"   "I"   "Q"   "V"   "V"   "V"  
    ## Mouse                         "N"   "I"   "Q"   "V"   "V"   "V"   "R"  
    ##                               [,24] [,25] [,26] [,27] [,28] [,29] [,30]
    ## Human                         "R"   "C"   "R"   "P"   "F"   "N"   "L"  
    ## Three-spined_stickleback_fish "E"   "R"   "G"   "S"   "A"   "K"   "R"  
    ## Clownfish                     "C"   "R"   "P"   "F"   "N"   "T"   "V"  
    ## Gilt-head_bream_fish          "C"   "R"   "P"   "F"   "N"   "T"   "M"  
    ## House_cat                     "R"   "C"   "R"   "P"   "F"   "N"   "L"  
    ## Mouse                         "C"   "R"   "P"   "F"   "N"   "L"   "A"  
    ##                               [,31] [,32] [,33] [,34] [,35] [,36] [,37]
    ## Human                         "A"   "E"   "R"   "K"   "A"   "S"   "A"  
    ## Three-spined_stickleback_fish "R"   "T"   "A"   "S"   "T"   "L"   "M"  
    ## Clownfish                     "E"   "R"   "K"   "S"   "S"   "Y"   "G"  
    ## Gilt-head_bream_fish          "E"   "R"   "K"   "S"   "S"   "Y"   "G"  
    ## House_cat                     "A"   "E"   "R"   "K"   "A"   "S"   "A"  
    ## Mouse                         "E"   "R"   "K"   "A"   "N"   "A"   "H"  
    ##                               [,38] [,39] [,40] [,41] [,42] [,43] [,44]
    ## Human                         "H"   "S"   "I"   "V"   "E"   "C"   "D"  
    ## Three-spined_stickleback_fish "N"   "A"   "Y"   "S"   "S"   "R"   "S"  
    ## Clownfish                     "V"   "I"   "D"   "C"   "D"   "Q"   "S"  
    ## Gilt-head_bream_fish          "V"   "I"   "D"   "C"   "D"   "S"   "N"  
    ## House_cat                     "H"   "S"   "V"   "V"   "E"   "C"   "D"  
    ## Mouse                         "S"   "V"   "V"   "E"   "C"   "D"   "H"  
    ##                               [,45] [,46] [,47] [,48] [,49] [,50] [,51]
    ## Human                         "P"   "V"   "R"   "K"   "E"   "V"   "S"  
    ## Three-spined_stickleback_fish "H"   "S"   "V"   "F"   "S"   "V"   "T"  
    ## Clownfish                     "R"   "K"   "E"   "V"   "M"   "V"   "K"  
    ## Gilt-head_bream_fish          "R"   "K"   "E"   "V"   "M"   "V"   "K"  
    ## House_cat                     "H"   "V"   "R"   "K"   "E"   "V"   "S"  
    ## Mouse                         "A"   "R"   "K"   "E"   "V"   "S"   "V"  
    ##                               [,52] [,53] [,54] [,55] [,56] [,57] [,58]
    ## Human                         "V"   "R"   "T"   "G"   "G"   "L"   "A"  
    ## Three-spined_stickleback_fish "I"   "H"   "M"   "K"   "E"   "I"   "T"  
    ## Clownfish                     "T"   "G"   "G"   "M"   "N"   "D"   "K"  
    ## Gilt-head_bream_fish          "T"   "G"   "G"   "M"   "N"   "D"   "K"  
    ## House_cat                     "V"   "R"   "T"   "G"   "G"   "L"   "A"  
    ## Mouse                         "R"   "T"   "A"   "G"   "L"   "T"   "D"  
    ##                               [,59] [,60] [,61] [,62] [,63] [,64] [,65]
    ## Human                         "D"   "K"   "S"   "S"   "R"   "K"   "T"  
    ## Three-spined_stickleback_fish "M"   "D"   "G"   "E"   "E"   "L"   "V"  
    ## Clownfish                     "A"   "S"   "R"   "K"   "T"   "Y"   "T"  
    ## Gilt-head_bream_fish          "A"   "S"   "R"   "K"   "T"   "Y"   "T"  
    ## House_cat                     "D"   "K"   "S"   "S"   "R"   "K"   "T"  
    ## Mouse                         "K"   "T"   "S"   "K"   "K"   "T"   "Y"  
    ##                               [,66] [,67] [,68] [,69] [,70] [,71] [,72]
    ## Human                         "Y"   "T"   "F"   "D"   "M"   "V"   "F"  
    ## Three-spined_stickleback_fish "K"   "I"   "G"   "K"   "L"   "N"   "L"  
    ## Clownfish                     "F"   "D"   "M"   "V"   "F"   "G"   "P"  
    ## Gilt-head_bream_fish          "F"   "D"   "M"   "V"   "F"   "G"   "P"  
    ## House_cat                     "Y"   "T"   "F"   "D"   "M"   "V"   "F"  
    ## Mouse                         "T"   "F"   "D"   "M"   "V"   "F"   "G"  
    ##                               [,73] [,74] [,75] [,76] [,77] [,78] [,79]
    ## Human                         "G"   "A"   "S"   "T"   "K"   "Q"   "I"  
    ## Three-spined_stickleback_fish "V"   "D"   "L"   "A"   "G"   "S"   "E"  
    ## Clownfish                     "A"   "A"   "K"   "Q"   "I"   "D"   "V"  
    ## Gilt-head_bream_fish          "A"   "A"   "K"   "Q"   "I"   "E"   "V"  
    ## House_cat                     "G"   "A"   "S"   "T"   "K"   "Q"   "I"  
    ## Mouse                         "A"   "S"   "T"   "K"   "Q"   "I"   "D"  
    ##                               [,80] [,81] [,82] [,83] [,84] [,85] [,86]
    ## Human                         "D"   "V"   "Y"   "R"   "S"   "V"   "V"  
    ## Three-spined_stickleback_fish "N"   "I"   "G"   "R"   "S"   "G"   "A"  
    ## Clownfish                     "Y"   "R"   "S"   "V"   "V"   "C"   "P"  
    ## Gilt-head_bream_fish          "Y"   "R"   "N"   "V"   "V"   "C"   "P"  
    ## House_cat                     "D"   "V"   "Y"   "R"   "S"   "V"   "V"  
    ## Mouse                         "V"   "Y"   "R"   "S"   "V"   "V"   "C"  
    ##                               [,87] [,88] [,89] [,90] [,91] [,92] [,93]
    ## Human                         "C"   "P"   "I"   "L"   "D"   "E"   "V"  
    ## Three-spined_stickleback_fish "V"   "D"   "K"   "R"   "A"   "R"   "E"  
    ## Clownfish                     "I"   "L"   "D"   "E"   "V"   "I"   "M"  
    ## Gilt-head_bream_fish          "I"   "L"   "D"   "E"   "V"   "I"   "M"  
    ## House_cat                     "C"   "P"   "I"   "L"   "D"   "E"   "V"  
    ## Mouse                         "P"   "I"   "L"   "D"   "E"   "V"   "I"  
    ##                               [,94] [,95] [,96] [,97] [,98] [,99] [,100]
    ## Human                         "I"   "M"   "G"   "Y"   "N"   "C"   "T"   
    ## Three-spined_stickleback_fish "A"   "G"   "N"   "I"   "N"   "Q"   "S"   
    ## Clownfish                     "G"   "Y"   "N"   "C"   "T"   "V"   "F"   
    ## Gilt-head_bream_fish          "G"   "Y"   "N"   "C"   "T"   "V"   "F"   
    ## House_cat                     "I"   "M"   "G"   "Y"   "N"   "C"   "T"   
    ## Mouse                         "M"   "G"   "Y"   "N"   "C"   "T"   "I"   
    ##                               [,101] [,102] [,103] [,104] [,105] [,106]
    ## Human                         "I"    "F"    "A"    "Y"    "G"    "Q"   
    ## Three-spined_stickleback_fish "L"    "L"    "T"    "L"    "G"    "R"   
    ## Clownfish                     "A"    "Y"    "G"    "Q"    "T"    "G"   
    ## Gilt-head_bream_fish          "A"    "Y"    "G"    "Q"    "T"    "G"   
    ## House_cat                     "I"    "F"    "A"    "Y"    "G"    "Q"   
    ## Mouse                         "F"    "A"    "Y"    "G"    "Q"    "T"   
    ##                               [,107] [,108] [,109] [,110] [,111] [,112]
    ## Human                         "T"    "G"    "T"    "G"    "K"    "T"   
    ## Three-spined_stickleback_fish "V"    "I"    "T"    "A"    "L"    "V"   
    ## Clownfish                     "T"    "G"    "K"    "T"    "F"    "T"   
    ## Gilt-head_bream_fish          "T"    "G"    "K"    "T"    "F"    "T"   
    ## House_cat                     "T"    "G"    "T"    "G"    "K"    "T"   
    ## Mouse                         "G"    "T"    "G"    "K"    "T"    "F"   
    ##                               [,113] [,114] [,115] [,116] [,117] [,118]
    ## Human                         "F"    "T"    "M"    "E"    "G"    "E"   
    ## Three-spined_stickleback_fish "E"    "K"    "R"    "P"    "H"    "V"   
    ## Clownfish                     "M"    "E"    "G"    "E"    "R"    "S"   
    ## Gilt-head_bream_fish          "M"    "E"    "G"    "E"    "R"    "S"   
    ## House_cat                     "F"    "T"    "M"    "E"    "G"    "E"   
    ## Mouse                         "T"    "M"    "E"    "G"    "E"    "R"   
    ##                               [,119] [,120] [,121] [,122] [,123] [,124]
    ## Human                         "R"    "S"    "P"    "N"    "E"    "E"   
    ## Three-spined_stickleback_fish "P"    "Y"    "R"    "E"    "S"    "K"   
    ## Clownfish                     "P"    "D"    "G"    "E"    "F"    "T"   
    ## Gilt-head_bream_fish          "P"    "N"    "E"    "Q"    "F"    "T"   
    ## House_cat                     "R"    "S"    "P"    "N"    "E"    "E"   
    ## Mouse                         "S"    "P"    "N"    "E"    "V"    "Y"   
    ##                               [,125] [,126] [,127] [,128] [,129] [,130]
    ## Human                         "Y"    "T"    "W"    "E"    "E"    "D"   
    ## Three-spined_stickleback_fish "L"    "T"    "R"    "I"    "L"    "Q"   
    ## Clownfish                     "W"    "E"    "E"    "D"    "P"    "L"   
    ## Gilt-head_bream_fish          "W"    "E"    "E"    "D"    "P"    "L"   
    ## House_cat                     "Y"    "T"    "W"    "E"    "E"    "D"   
    ## Mouse                         "T"    "W"    "E"    "E"    "D"    "P"   
    ##                               [,131] [,132] [,133] [,134] [,135] [,136]
    ## Human                         "P"    "L"    "A"    "G"    "I"    "I"   
    ## Three-spined_stickleback_fish "D"    "S"    "L"    "G"    "G"    "R"   
    ## Clownfish                     "A"    "G"    "I"    "I"    "P"    "R"   
    ## Gilt-head_bream_fish          "A"    "G"    "I"    "I"    "P"    "R"   
    ## House_cat                     "P"    "L"    "A"    "G"    "I"    "I"   
    ## Mouse                         "L"    "A"    "G"    "I"    "I"    "P"   
    ##                               [,137] [,138] [,139] [,140] [,141] [,142]
    ## Human                         "P"    "R"    "T"    "L"    "H"    "Q"   
    ## Three-spined_stickleback_fish "T"    "K"    "T"    "S"    "I"    "I"   
    ## Clownfish                     "T"    "L"    "H"    "Q"    "I"    "F"   
    ## Gilt-head_bream_fish          "T"    "L"    "H"    "Q"    "I"    "F"   
    ## House_cat                     "P"    "R"    "T"    "L"    "H"    "Q"   
    ## Mouse                         "R"    "T"    "L"    "H"    "Q"    "I"   
    ##                               [,143] [,144] [,145] [,146] [,147] [,148]
    ## Human                         "I"    "F"    "E"    "K"    "L"    "T"   
    ## Three-spined_stickleback_fish "A"    "T"    "V"    "S"    "P"    "S"   
    ## Clownfish                     "E"    "K"    "L"    "S"    "E"    "N"   
    ## Gilt-head_bream_fish          "E"    "K"    "L"    "S"    "E"    "N"   
    ## House_cat                     "I"    "F"    "E"    "K"    "L"    "T"   
    ## Mouse                         "F"    "E"    "K"    "L"    "T"    "D"   
    ##                               [,149] [,150] [,151] [,152] [,153] [,154]
    ## Human                         "D"    "N"    "G"    "T"    "E"    "F"   
    ## Three-spined_stickleback_fish "S"    "S"    "N"    "L"    "E"    "E"   
    ## Clownfish                     "G"    "T"    "E"    "F"    "S"    "V"   
    ## Gilt-head_bream_fish          "G"    "T"    "E"    "F"    "S"    "V"   
    ## House_cat                     "D"    "N"    "G"    "T"    "E"    "F"   
    ## Mouse                         "N"    "G"    "T"    "E"    "F"    "S"   
    ##                               [,155] [,156] [,157] [,158] [,159] [,160]
    ## Human                         "S"    "V"    "K"    "V"    "S"    "L"   
    ## Three-spined_stickleback_fish "T"    "L"    "S"    "T"    "L"    "E"   
    ## Clownfish                     "K"    "V"    "S"    "L"    "L"    "E"   
    ## Gilt-head_bream_fish          "K"    "V"    "S"    "L"    "L"    "E"   
    ## House_cat                     "S"    "V"    "K"    "V"    "S"    "L"   
    ## Mouse                         "V"    "K"    "V"    "S"    "L"    "L"   
    ##                               [,161] [,162] [,163] [,164] [,165] [,166]
    ## Human                         "L"    "E"    "I"    "Y"    "N"    "E"   
    ## Three-spined_stickleback_fish "Y"    "A"    "S"    "R"    "A"    "K"   
    ## Clownfish                     "I"    "Y"    "N"    "E"    "E"    "L"   
    ## Gilt-head_bream_fish          "I"    "Y"    "N"    "E"    "E"    "L"   
    ## House_cat                     "L"    "E"    "I"    "Y"    "N"    "E"   
    ## Mouse                         "E"    "I"    "Y"    "N"    "E"    "E"   
    ##                               [,167] [,168] [,169] [,170] [,171] [,172]
    ## Human                         "E"    "L"    "F"    "D"    "L"    "L"   
    ## Three-spined_stickleback_fish "N"    "I"    "M"    "N"    "K"    "P"   
    ## Clownfish                     "F"    "D"    "L"    "L"    "S"    "P"   
    ## Gilt-head_bream_fish          "F"    "D"    "L"    "L"    "S"    "P"   
    ## House_cat                     "E"    "L"    "F"    "D"    "L"    "L"   
    ## Mouse                         "L"    "F"    "D"    "L"    "L"    "S"   
    ##                               [,173] [,174] [,175] [,176] [,177] [,178]
    ## Human                         "N"    "P"    "S"    "S"    "D"    "V"   
    ## Three-spined_stickleback_fish "E"    "V"    "N"    "Q"    "K"    "L"   
    ## Clownfish                     "T"    "E"    "D"    "V"    "N"    "E"   
    ## Gilt-head_bream_fish          "T"    "E"    "D"    "V"    "N"    "E"   
    ## House_cat                     "N"    "P"    "S"    "S"    "D"    "V"   
    ## Mouse                         "P"    "S"    "S"    "D"    "V"    "S"   
    ##                               [,179] [,180] [,181] [,182] [,183] [,184]
    ## Human                         "S"    "E"    "R"    "L"    "Q"    "M"   
    ## Three-spined_stickleback_fish "T"    "K"    "R"    "T"    "L"    "I"   
    ## Clownfish                     "R"    "L"    "Q"    "L"    "F"    "D"   
    ## Gilt-head_bream_fish          "R"    "L"    "Q"    "L"    "F"    "D"   
    ## House_cat                     "S"    "E"    "R"    "L"    "Q"    "M"   
    ## Mouse                         "E"    "R"    "L"    "Q"    "M"    "F"   
    ##                               [,185] [,186] [,187] [,188] [,189] [,190]
    ## Human                         "F"    "D"    "D"    "P"    "R"    "N"   
    ## Three-spined_stickleback_fish "K"    "E"    "Y"    "T"    "E"    "E"   
    ## Clownfish                     "D"    "P"    "R"    "N"    "K"    "R"   
    ## Gilt-head_bream_fish          "D"    "P"    "R"    "N"    "K"    "R"   
    ## House_cat                     "F"    "D"    "D"    "P"    "R"    "N"   
    ## Mouse                         "D"    "D"    "P"    "R"    "N"    "K"   
    ##                               [,191] [,192] [,193] [,194] [,195] [,196]
    ## Human                         "K"    "R"    "G"    "V"    "I"    "I"   
    ## Three-spined_stickleback_fish "I"    "E"    "R"    "L"    "K"    "R"   
    ## Clownfish                     "G"    "V"    "V"    "V"    "K"    "G"   
    ## Gilt-head_bream_fish          "G"    "V"    "V"    "V"    "K"    "G"   
    ## House_cat                     "K"    "R"    "G"    "V"    "I"    "I"   
    ## Mouse                         "R"    "G"    "V"    "I"    "I"    "K"   
    ##                               [,197] [,198] [,199] [,200] [,201] [,202]
    ## Human                         "K"    "G"    "L"    "E"    "E"    "I"   
    ## Three-spined_stickleback_fish "D"    "L"    "A"    "A"    "T"    "R"   
    ## Clownfish                     "L"    "E"    "E"    "V"    "T"    "V"   
    ## Gilt-head_bream_fish          "L"    "E"    "E"    "V"    "T"    "V"   
    ## House_cat                     "K"    "G"    "L"    "E"    "E"    "I"   
    ## Mouse                         "G"    "L"    "E"    "E"    "I"    "T"   
    ##                               [,203] [,204] [,205] [,206] [,207] [,208]
    ## Human                         "T"    "V"    "H"    "N"    "K"    "D"   
    ## Three-spined_stickleback_fish "D"    "K"    "N"    "G"    "I"    "Y"   
    ## Clownfish                     "H"    "N"    "K"    "D"    "E"    "V"   
    ## Gilt-head_bream_fish          "H"    "N"    "K"    "D"    "E"    "V"   
    ## House_cat                     "T"    "V"    "H"    "N"    "K"    "D"   
    ## Mouse                         "V"    "H"    "N"    "K"    "D"    "E"   
    ##                               [,209] [,210] [,211] [,212] [,213] [,214]
    ## Human                         "E"    "V"    "Y"    "Q"    "I"    "L"   
    ## Three-spined_stickleback_fish "L"    "S"    "A"    "E"    "N"    "Y"   
    ## Clownfish                     "Y"    "Q"    "I"    "L"    "E"    "R"   
    ## Gilt-head_bream_fish          "Y"    "Q"    "I"    "L"    "E"    "R"   
    ## House_cat                     "E"    "V"    "Y"    "Q"    "I"    "L"   
    ## Mouse                         "V"    "Y"    "Q"    "I"    "L"    "E"   
    ##                               [,215] [,216] [,217] [,218] [,219] [,220]
    ## Human                         "E"    "K"    "G"    "A"    "A"    "K"   
    ## Three-spined_stickleback_fish "E"    "S"    "M"    "M"    "G"    "Q"   
    ## Clownfish                     "G"    "A"    "A"    "K"    "R"    "R"   
    ## Gilt-head_bream_fish          "G"    "S"    "A"    "K"    "R"    "R"   
    ## House_cat                     "E"    "K"    "G"    "A"    "A"    "K"   
    ## Mouse                         "K"    "G"    "A"    "A"    "K"    "R"   
    ##                               [,221] [,222] [,223] [,224] [,225] [,226]
    ## Human                         "R"    "T"    "T"    "A"    "A"    "T"   
    ## Three-spined_stickleback_fish "I"    "T"    "S"    "H"    "E"    "V"   
    ## Clownfish                     "T"    "A"    "S"    "T"    "L"    "M"   
    ## Gilt-head_bream_fish          "T"    "A"    "S"    "T"    "L"    "M"   
    ## House_cat                     "R"    "T"    "T"    "A"    "A"    "T"   
    ## Mouse                         "T"    "T"    "A"    "A"    "T"    "L"   
    ##                               [,227] [,228] [,229] [,230] [,231] [,232]
    ## Human                         "L"    "M"    "N"    "A"    "Y"    "S"   
    ## Three-spined_stickleback_fish "H"    "T"    "V"    "E"    "Y"    "S"   
    ## Clownfish                     "N"    "A"    "Y"    "S"    "S"    "R"   
    ## Gilt-head_bream_fish          "N"    "A"    "Y"    "S"    "S"    "R"   
    ## House_cat                     "L"    "M"    "N"    "A"    "Y"    "S"   
    ## Mouse                         "M"    "N"    "A"    "Y"    "S"    "S"   
    ##                               [,233] [,234] [,235] [,236] [,237] [,238]
    ## Human                         "S"    "R"    "S"    "H"    "S"    "V"   
    ## Three-spined_stickleback_fish "D"    "R"    "I"    "A"    "A"    "M"   
    ## Clownfish                     "S"    "H"    "S"    "V"    "F"    "S"   
    ## Gilt-head_bream_fish          "S"    "H"    "S"    "V"    "F"    "S"   
    ## House_cat                     "S"    "R"    "S"    "H"    "S"    "V"   
    ## Mouse                         "R"    "S"    "H"    "S"    "V"    "F"   
    ##                               [,239] [,240] [,241] [,242] [,243] [,244]
    ## Human                         "F"    "S"    "V"    "T"    "I"    "H"   
    ## Three-spined_stickleback_fish "E"    "E"    "E"    "I"    "K"    "K"   
    ## Clownfish                     "V"    "T"    "I"    "H"    "M"    "K"   
    ## Gilt-head_bream_fish          "V"    "T"    "I"    "H"    "M"    "K"   
    ## House_cat                     "F"    "S"    "V"    "T"    "I"    "H"   
    ## Mouse                         "S"    "V"    "T"    "I"    "H"    "M"   
    ##                               [,245] [,246] [,247] [,248] [,249] [,250]
    ## Human                         "M"    "K"    "E"    "T"    "T"    "I"   
    ## Three-spined_stickleback_fish "V"    "T"    "E"    "L"    "F"    "V"   
    ## Clownfish                     "E"    "I"    "T"    "V"    "D"    "G"   
    ## Gilt-head_bream_fish          "E"    "I"    "T"    "L"    "D"    "G"   
    ## House_cat                     "M"    "K"    "E"    "T"    "T"    "I"   
    ## Mouse                         "K"    "E"    "T"    "T"    "I"    "D"   
    ##                               [,251] [,252] [,253] [,254] [,255] [,256]
    ## Human                         "D"    "G"    "E"    "E"    "L"    "V"   
    ## Three-spined_stickleback_fish "D"    "S"    "K"    "T"    "R"    "L"   
    ## Clownfish                     "E"    "E"    "L"    "V"    "K"    "I"   
    ## Gilt-head_bream_fish          "E"    "E"    "L"    "V"    "K"    "I"   
    ## House_cat                     "D"    "G"    "E"    "E"    "L"    "V"   
    ## Mouse                         "G"    "E"    "E"    "L"    "V"    "K"   
    ##                               [,257] [,258] [,259] [,260] [,261] [,262]
    ## Human                         "K"    "I"    "G"    "K"    "L"    "N"   
    ## Three-spined_stickleback_fish "E"    "L"    "C"    "A"    "V"    "D"   
    ## Clownfish                     "G"    "K"    "L"    "N"    "L"    "V"   
    ## Gilt-head_bream_fish          "G"    "K"    "L"    "N"    "L"    "V"   
    ## House_cat                     "K"    "I"    "G"    "K"    "L"    "N"   
    ## Mouse                         "I"    "G"    "K"    "L"    "N"    "L"   
    ##                               [,263] [,264] [,265] [,266] [,267] [,268]
    ## Human                         "L"    "V"    "D"    "L"    "A"    "G"   
    ## Three-spined_stickleback_fish "L"    "D"    "E"    "K"    "Q"    "Q"   
    ## Clownfish                     "D"    "L"    "A"    "G"    "S"    "E"   
    ## Gilt-head_bream_fish          "D"    "L"    "A"    "G"    "S"    "E"   
    ## House_cat                     "L"    "V"    "D"    "L"    "A"    "G"   
    ## Mouse                         "V"    "D"    "L"    "A"    "G"    "S"   
    ##                               [,269] [,270] [,271] [,272] [,273] [,274]
    ## Human                         "S"    "E"    "N"    "I"    "G"    "R"   
    ## Three-spined_stickleback_fish "R"    "L"    "E"    "E"    "T"    "S"   
    ## Clownfish                     "N"    "I"    "G"    "R"    "S"    "G"   
    ## Gilt-head_bream_fish          "N"    "I"    "G"    "R"    "S"    "G"   
    ## House_cat                     "S"    "E"    "N"    "I"    "G"    "R"   
    ## Mouse                         "E"    "N"    "I"    "G"    "R"    "S"   
    ##                               [,275] [,276] [,277] [,278] [,279] [,280]
    ## Human                         "S"    "G"    "A"    "V"    "D"    "K"   
    ## Three-spined_stickleback_fish "R"    "D"    "L"    "Q"    "H"    "T"   
    ## Clownfish                     "A"    "V"    "D"    "K"    "R"    "A"   
    ## Gilt-head_bream_fish          "A"    "V"    "D"    "K"    "R"    "A"   
    ## House_cat                     "S"    "G"    "A"    "V"    "D"    "K"   
    ## Mouse                         "G"    "A"    "V"    "D"    "K"    "R"   
    ##                               [,281] [,282] [,283] [,284] [,285] [,286]
    ## Human                         "R"    "A"    "R"    "E"    "A"    "G"   
    ## Three-spined_stickleback_fish "K"    "E"    "K"    "L"    "M"    "E"   
    ## Clownfish                     "R"    "E"    "A"    "G"    "N"    "I"   
    ## Gilt-head_bream_fish          "R"    "E"    "A"    "G"    "N"    "I"   
    ## House_cat                     "R"    "A"    "R"    "E"    "A"    "G"   
    ## Mouse                         "A"    "R"    "E"    "A"    "G"    "N"   
    ##                               [,287] [,288] [,289] [,290] [,291] [,292]
    ## Human                         "N"    "I"    "N"    "Q"    "S"    "L"   
    ## Three-spined_stickleback_fish "X"    "E"    "F"    "V"    "C"    "S"   
    ## Clownfish                     "N"    "Q"    "S"    "L"    "L"    "T"   
    ## Gilt-head_bream_fish          "N"    "Q"    "S"    "L"    "L"    "T"   
    ## House_cat                     "N"    "I"    "N"    "Q"    "S"    "L"   
    ## Mouse                         "I"    "N"    "Q"    "S"    "L"    "L"   
    ##                               [,293] [,294] [,295] [,296] [,297] [,298]
    ## Human                         "L"    "T"    "L"    "G"    "R"    "V"   
    ## Three-spined_stickleback_fish "E"    "L"    "T"    "L"    "V"    "Q"   
    ## Clownfish                     "L"    "G"    "R"    "V"    "I"    "T"   
    ## Gilt-head_bream_fish          "L"    "G"    "R"    "V"    "I"    "T"   
    ## House_cat                     "L"    "T"    "L"    "G"    "R"    "V"   
    ## Mouse                         "T"    "L"    "G"    "R"    "V"    "I"   
    ##                               [,299] [,300] [,301] [,302] [,303] [,304]
    ## Human                         "I"    "T"    "A"    "L"    "V"    "E"   
    ## Three-spined_stickleback_fish "E"    "S"    "L"    "Y"    "D"    "T"   
    ## Clownfish                     "A"    "L"    "V"    "E"    "K"    "R"   
    ## Gilt-head_bream_fish          "A"    "L"    "V"    "E"    "K"    "R"   
    ## House_cat                     "I"    "T"    "A"    "L"    "V"    "E"   
    ## Mouse                         "T"    "A"    "L"    "V"    "E"    "R"   
    ##                               [,305] [,306] [,307] [,308] [,309] [,310]
    ## Human                         "R"    "T"    "P"    "H"    "V"    "P"   
    ## Three-spined_stickleback_fish "A"    "G"    "R"    "L"    "L"    "S"   
    ## Clownfish                     "P"    "H"    "I"    "P"    "Y"    "R"   
    ## Gilt-head_bream_fish          "P"    "H"    "V"    "P"    "Y"    "R"   
    ## House_cat                     "R"    "T"    "P"    "H"    "V"    "P"   
    ## Mouse                         "T"    "P"    "H"    "I"    "P"    "Y"   
    ##                               [,311] [,312] [,313] [,314] [,315] [,316]
    ## Human                         "Y"    "R"    "E"    "S"    "K"    "L"   
    ## Three-spined_stickleback_fish "T"    "V"    "D"    "A"    "S"    "T"   
    ## Clownfish                     "E"    "S"    "K"    "L"    "T"    "R"   
    ## Gilt-head_bream_fish          "E"    "S"    "K"    "L"    "T"    "R"   
    ## House_cat                     "Y"    "R"    "E"    "S"    "K"    "L"   
    ## Mouse                         "R"    "E"    "S"    "K"    "L"    "T"   
    ##                               [,317] [,318] [,319] [,320] [,321] [,322]
    ## Human                         "T"    "R"    "I"    "L"    "Q"    "D"   
    ## Three-spined_stickleback_fish "G"    "D"    "V"    "C"    "G"    "L"   
    ## Clownfish                     "I"    "L"    "Q"    "D"    "S"    "L"   
    ## Gilt-head_bream_fish          "I"    "L"    "Q"    "D"    "S"    "L"   
    ## House_cat                     "T"    "R"    "I"    "L"    "Q"    "D"   
    ## Mouse                         "R"    "I"    "L"    "Q"    "D"    "S"   
    ##                               [,323] [,324] [,325] [,326] [,327] [,328]
    ## Human                         "S"    "L"    "G"    "G"    "R"    "T"   
    ## Three-spined_stickleback_fish "P"    "G"    "Q"    "L"    "D"    "R"   
    ## Clownfish                     "G"    "G"    "R"    "T"    "K"    "T"   
    ## Gilt-head_bream_fish          "G"    "G"    "R"    "T"    "K"    "T"   
    ## House_cat                     "S"    "L"    "G"    "G"    "R"    "T"   
    ## Mouse                         "L"    "G"    "G"    "R"    "T"    "R"   
    ##                               [,329] [,330] [,331] [,332] [,333] [,334]
    ## Human                         "R"    "T"    "S"    "I"    "I"    "A"   
    ## Three-spined_stickleback_fish "X"    "K"    "X"    "V"    "E"    "Q"   
    ## Clownfish                     "S"    "I"    "I"    "A"    "T"    "V"   
    ## Gilt-head_bream_fish          "S"    "I"    "I"    "A"    "T"    "V"   
    ## House_cat                     "R"    "T"    "S"    "I"    "I"    "A"   
    ## Mouse                         "T"    "S"    "I"    "I"    "A"    "T"   
    ##                               [,335] [,336] [,337] [,338] [,339] [,340]
    ## Human                         "T"    "I"    "S"    "P"    "A"    "S"   
    ## Three-spined_stickleback_fish "H"    "Y"    "S"    "G"    "V"    "Q"   
    ## Clownfish                     "S"    "P"    "S"    "S"    "S"    "N"   
    ## Gilt-head_bream_fish          "S"    "P"    "S"    "S"    "S"    "N"   
    ## House_cat                     "T"    "I"    "S"    "P"    "A"    "S"   
    ## Mouse                         "I"    "S"    "P"    "A"    "S"    "F"   
    ##                               [,341] [,342] [,343] [,344] [,345] [,346]
    ## Human                         "L"    "N"    "L"    "E"    "E"    "T"   
    ## Three-spined_stickleback_fish "Q"    "S"    "S"    "L"    "S"    "A"   
    ## Clownfish                     "L"    "E"    "E"    "T"    "L"    "S"   
    ## Gilt-head_bream_fish          "L"    "E"    "E"    "T"    "L"    "S"   
    ## House_cat                     "L"    "N"    "L"    "E"    "E"    "T"   
    ## Mouse                         "N"    "L"    "E"    "E"    "T"    "L"   
    ##                               [,347] [,348]
    ## Human                         "L"    "S"   
    ## Three-spined_stickleback_fish "W"    "X"   
    ## Clownfish                     "T"    "L"   
    ## Gilt-head_bream_fish          "T"    "L"   
    ## House_cat                     "L"    "S"   
    ## Mouse                         "S"    "T"

``` r
seqidentity(x, normalize=TRUE, similarity=FALSE, ncore=1, nseg.scale=1)
```

    ##                               Human Three-spined_stickleback_fish
    ## Human                         1.000                         0.057
    ## Three-spined_stickleback_fish 0.057                         1.000
    ## Clownfish                     0.063                         0.069
    ## Gilt-head_bream_fish          0.066                         0.069
    ## House_cat                     0.895                         0.060
    ## Mouse                         0.080                         0.055
    ## Yellow_perch_fish             0.070                         0.069
    ## Zander_fish                   0.071                         0.072
    ## Zebrafish                     0.063                         0.057
    ##                               Clownfish Gilt-head_bream_fish House_cat
    ## Human                             0.063                0.066     0.895
    ## Three-spined_stickleback_fish     0.069                0.069     0.060
    ## Clownfish                         1.000                0.835     0.075
    ## Gilt-head_bream_fish              0.835                1.000     0.079
    ## House_cat                         0.075                0.079     1.000
    ## Mouse                             0.099                0.091     0.090
    ## Yellow_perch_fish                 0.807                0.802     0.082
    ## Zander_fish                       0.816                0.803     0.083
    ## Zebrafish                         0.120                0.113     0.071
    ##                               Mouse Yellow_perch_fish Zander_fish
    ## Human                         0.080             0.070       0.071
    ## Three-spined_stickleback_fish 0.055             0.069       0.072
    ## Clownfish                     0.099             0.807       0.816
    ## Gilt-head_bream_fish          0.091             0.802       0.803
    ## House_cat                     0.090             0.082       0.083
    ## Mouse                         1.000             0.091       0.093
    ## Yellow_perch_fish             0.091             1.000       0.907
    ## Zander_fish                   0.093             0.907       1.000
    ## Zebrafish                     0.451             0.095       0.097
    ##                               Zebrafish
    ## Human                             0.063
    ## Three-spined_stickleback_fish     0.057
    ## Clownfish                         0.120
    ## Gilt-head_bream_fish              0.113
    ## House_cat                         0.071
    ## Mouse                             0.451
    ## Yellow_perch_fish                 0.095
    ## Zander_fish                       0.097
    ## Zebrafish                         1.000

``` r
y<- seqidentity(x, normalize=TRUE, similarity=FALSE, ncore=1, nseg.scale=1)
heatmap(y, margins = c(10,10))
```

![](FAGP_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
png( "~/Desktop/analines_R_plot_saved_to_png_december_5_2019_bioinformatics_fall_2020_final_project_verion_1.1.1.0.png")
heatmap(y, margins = c(12,12), cexRow = 1, cexCol = 1)
z <- heatmap(y, margins = c(12,12), cexRow = 1, cexCol = 1)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
rowSums(y)
```

    ##                         Human Three-spined_stickleback_fish 
    ##                         2.365                         1.508 
    ##                     Clownfish          Gilt-head_bream_fish 
    ##                         3.884                         3.858 
    ##                     House_cat                         Mouse 
    ##                         2.435                         2.050 
    ##             Yellow_perch_fish                   Zander_fish 
    ##                         3.923                         3.942 
    ##                     Zebrafish 
    ##                         2.067

``` r
# zander fish is bomb
```

``` r
fish <- read.fasta("Zander_fish copy.txt")
results <- blast.pdb(fish)
```

    ##  Searching ... please wait (updates every 5 seconds) RID = YWWTRSBY014 
    ##  ...
    ##  Reporting 101 hits

``` r
results
```

    ## $hit.tbl
    ##         queryid subjectids identity alignmentlength mismatches gapopens
    ## 1   Query_47107     3WPN_A   85.366             369         52        2
    ## 2   Query_47107     3HQD_A   85.366             369         52        2
    ## 3   Query_47107     1II6_A   85.326             368         52        2
    ## 4   Query_47107     4ZCA_A   85.095             369         53        2
    ## 5   Query_47107     4ZHI_A   85.095             369         53        2
    ## 6   Query_47107     4AP0_A   85.326             368         52        2
    ## 7   Query_47107     1Q0B_A   85.675             363         51        1
    ## 8   Query_47107     4A1Z_A   85.054             368         53        2
    ## 9   Query_47107     4A28_A   85.054             368         53        2
    ## 10  Query_47107     1X88_A   86.313             358         48        1
    ## 11  Query_47107     4B7B_A   84.783             368         54        2
    ## 12  Query_47107     5ZO7_A   86.686             353         46        1
    ## 13  Query_47107     4AQV_C   83.924             367         57        2
    ## 14  Query_47107     3ZCW_A   86.494             348         46        1
    ## 15  Query_47107     2WBE_C   56.383             376        152        6
    ## 16  Query_47107     5MM4_K   52.332             386        145        8
    ## 17  Query_47107     5MM7_K   52.030             394        149        9
    ## 18  Query_47107     5M5I_C   50.000             368        165        6
    ## 19  Query_47107     5MLV_G   50.000             368        165        6
    ## 20  Query_47107     6S8M_K   50.000             368        165        6
    ## 21  Query_47107     3B6U_A   46.832             363        175        8
    ## 22  Query_47107     3B6V_A   41.538             390        180       10
    ## 23  Query_47107     1GOJ_A   44.321             361        174        8
    ## 24  Query_47107     2VVG_A   42.297             357        180        8
    ## 25  Query_47107     3ZFC_A   43.094             362        182        9
    ## 26  Query_47107     3ZFD_A   43.094             362        182        9
    ## 27  Query_47107     6MLQ_C   41.274             361        182       10
    ## 28  Query_47107     4BN2_A   44.384             365        172       12
    ## 29  Query_47107     2Y5W_A   41.711             374        187       10
    ## 30  Query_47107     6A20_A   38.520             392        204       10
    ## 31  Query_47107     6A1Z_A   38.520             392        204       10
    ## 32  Query_47107     4HNA_K   42.254             355        176        8
    ## 33  Query_47107     3J8X_K   42.254             355        176        8
    ## 34  Query_47107     4A14_A   41.143             350        176       10
    ## 35  Query_47107     1MKJ_A   42.254             355        176        8
    ## 36  Query_47107     4ATX_C   42.254             355        176        8
    ## 37  Query_47107     4UXT_C   42.017             357        173       11
    ## 38  Query_47107     3LRE_A   39.669             363        181        9
    ## 39  Query_47107     3J6H_K   41.972             355        177        9
    ## 40  Query_47107     5ZBS_B   38.462             377        196        9
    ## 41  Query_47107     4LNU_K   42.151             344        170        8
    ## 42  Query_47107     5OAM_K   39.669             363        181        9
    ## 43  Query_47107     5ZBR_B   38.462             377        196        9
    ## 44  Query_47107     1BG2_A   42.151             344        170        8
    ## 45  Query_47107     2XT3_A   40.571             350        178       10
    ## 46  Query_47107     5LT3_A   41.860             344        171        8
    ## 47  Query_47107     5LT1_A   41.860             344        171        8
    ## 48  Query_47107     3WRD_A   41.880             351        175        9
    ## 49  Query_47107     5HNY_K   41.737             357        174       11
    ## 50  Query_47107     1T5C_A   39.943             353        185        9
    ## 51  Query_47107     5HNW_K   41.176             357        176       11
    ## 52  Query_47107     1SDM_A   37.726             387        205       12
    ## 53  Query_47107     5GSZ_A   37.705             366        197        9
    ## 54  Query_47107     5GSY_K   37.705             366        197        9
    ## 55  Query_47107     3H4S_A   38.920             352        191       10
    ## 56  Query_47107     3GBJ_A   38.333             360        186        9
    ## 57  Query_47107     1IA0_K   36.220             381        198       10
    ## 58  Query_47107     1I6I_A   35.979             378        197       10
    ## 59  Query_47107     1VFV_A   35.979             378        197       10
    ## 60  Query_47107     4OZQ_A   40.220             363        178       11
    ## 61  Query_47107     1I5S_A   36.074             377        196       10
    ## 62  Query_47107     4UXO_C   35.544             377        198       10
    ## 63  Query_47107     5WDE_A   38.857             350        179        9
    ## 64  Query_47107     3NWN_A   38.961             308        176        5
    ## 65  Query_47107     6NJE_A   35.556             360        207        8
    ## 66  Query_47107     2OWM_A   39.650             343        151       11
    ## 67  Query_47107     2NCD_A   35.440             364        201       12
    ## 68  Query_47107     1CZ7_A   35.440             364        201       12
    ## 69  Query_47107     3U06_A   35.376             359        199       11
    ## 70  Query_47107     5W3D_A   35.440             364        201       12
    ## 71  Query_47107     3L1C_A   35.632             348        191       11
    ## 72  Query_47107     4ETP_A   36.188             362        180       10
    ## 73  Query_47107     1N6M_A   35.165             364        202       12
    ## 74  Query_47107     3PXN_A   34.688             320        180        7
    ## 75  Query_47107     3DC4_A   34.906             318        178        7
    ## 76  Query_47107     1F9T_A   35.027             374        168       10
    ## 77  Query_47107     3KAR_A   35.027             374        168       10
    ## 78  Query_47107     1F9W_A   34.759             374        169       10
    ## 79  Query_47107     4GKR_A   36.570             309        160        8
    ## 80  Query_47107     1F9U_A   34.759             374        169       10
    ## 81  Query_47107     1F9V_A   34.759             374        169       10
    ## 82  Query_47107     4Y05_A   36.494             348        189       11
    ## 83  Query_47107     2HEH_A   35.447             347        200        9
    ## 84  Query_47107     3EDL_D   34.870             347        202        9
    ## 85  Query_47107     3T0Q_A   36.415             357        187       13
    ## 86  Query_47107     5MIO_C   35.795             352        192       11
    ## 87  Query_47107     1V8J_A   34.870             347        202        9
    ## 88  Query_47107     5XJA_A   36.242             298        172        6
    ## 89  Query_47107     5WDH_A   36.585             369        183       14
    ## 90  Query_47107     2GRY_A   34.393             346        205        7
    ## 91  Query_47107     4H1G_A   35.294             323        183        8
    ## 92  Query_47107     6B0I_K   36.513             304        167        9
    ## 93  Query_47107     3J2U_K   36.513             304        167        9
    ## 94  Query_47107     1RY6_A   33.903             351        198       11
    ## 95  Query_47107     5ND2_C   40.625             224        114        7
    ## 96  Query_47107     5ND2_C   26.087             138         81        4
    ## 97  Query_47107     2KIN_A   37.743             257        132        8
    ## 98  Query_47107     5X3E_A   29.551             379        199       11
    ## 99  Query_47107     2KIN_B   55.172              87         38        1
    ## 100 Query_47107     3KIN_B   55.422              83         36        1
    ## 101 Query_47107     5JV3_A   71.429              35         10        0
    ##     q.start q.end s.start s.end    evalue bitscore positives mlog.evalue
    ## 1         1   367       7   375  0.00e+00    677.0     94.58   709.19621
    ## 2         1   367       1   369  0.00e+00    676.0     94.58   709.19621
    ## 3         1   366       1   368  0.00e+00    675.0     94.57   709.19621
    ## 4         1   367       1   369  0.00e+00    675.0     94.58   709.19621
    ## 5         1   367       1   369  0.00e+00    674.0     94.31   709.19621
    ## 6         1   366       3   370  0.00e+00    674.0     94.57   709.19621
    ## 7         5   366       5   367  0.00e+00    674.0     95.04   709.19621
    ## 8         1   366       1   368  0.00e+00    672.0     94.29   709.19621
    ## 9         1   366       1   368  0.00e+00    672.0     94.29   709.19621
    ## 10       10   366       2   359  0.00e+00    671.0     95.81   709.19621
    ## 11        1   366       1   368  0.00e+00    668.0     94.02   709.19621
    ## 12       16   367       8   360  0.00e+00    664.0     95.75   709.19621
    ## 13        1   365       1   367  0.00e+00    658.0     93.19   709.19621
    ## 14       15   361       1   348  0.00e+00    654.0     95.69   709.19621
    ## 15        1   369       3   373 3.21e-136    416.0     74.73   311.98530
    ## 16       17   368       5   385 4.11e-113    356.0     67.88   258.77869
    ## 17       10   368      69   457 4.61e-113    358.0     67.51   258.66389
    ## 18       12   365       4   366 1.25e-106    338.0     67.93   243.85088
    ## 19       12   365      67   429 2.39e-105    337.0     67.93   240.90014
    ## 20       12   365      73   435 6.10e-105    336.0     67.93   239.96315
    ## 21       18   376      23   371  4.11e-95    307.0     63.91   217.33216
    ## 22       18   376      23   395  1.13e-84    280.0     59.23   193.29493
    ## 23       17   372       7   345  7.03e-82    271.0     61.22   186.86179
    ## 24       17   365       5   343  3.60e-81    269.0     61.06   185.22846
    ## 25        9   364       1   344  2.67e-79    264.0     61.33   180.92214
    ## 26        9   364       1   344  2.76e-79    264.0     61.33   180.88899
    ## 27       18   366      17   359  1.35e-69    239.0     59.56   158.57827
    ## 28       15   365      10   357  5.33e-69    236.0     56.71   157.20502
    ## 29       18   385      13   361  3.40e-68    234.0     57.49   155.35201
    ## 30       18   387       6   382  4.49e-68    236.0     56.38   155.07393
    ## 31       18   387       6   382  1.16e-67    233.0     56.38   154.12478
    ## 32       17   368       8   336  3.99e-67    230.0     55.21   152.88941
    ## 33       17   368       8   336  8.63e-67    229.0     55.21   152.11796
    ## 34       18   355      13   344  1.52e-66    228.0     59.43   151.55191
    ## 35       17   368       8   336  2.64e-66    228.0     55.21   150.99984
    ## 36       17   368       8   336  4.99e-66    227.0     54.93   150.36318
    ## 37       18   368      10   338  1.67e-65    225.0     56.30   149.15521
    ## 38       16   357      10   355  5.79e-65    224.0     56.75   147.91190
    ## 39       18   368       9   338  6.04e-65    224.0     55.77   147.86963
    ## 40       18   372       4   366  6.36e-65    225.0     55.70   147.81800
    ## 41       17   357       8   325  8.37e-65    223.0     55.23   147.54338
    ## 42       16   357      13   358  1.10e-64    224.0     56.75   147.27014
    ## 43       18   372       4   366  1.29e-64    224.0     55.70   147.11080
    ## 44       17   357       8   325  2.93e-64    221.0     55.23   146.29044
    ## 45       18   355      13   344  3.34e-64    222.0     59.14   146.15948
    ## 46       17   357       8   325  4.57e-64    221.0     54.94   145.84593
    ## 47       17   357       8   325  5.25e-64    220.0     54.94   145.70722
    ## 48       18   364       9   334  7.15e-64    221.0     55.56   145.39833
    ## 49       17   367      24   352  1.10e-63    221.0     56.30   144.96755
    ## 50       18   365       6   336  4.56e-63    219.0     55.52   143.54554
    ## 51       17   367      24   352  2.07e-62    218.0     56.02   142.03273
    ## 52       12   395       2   355  3.20e-62    217.0     56.59   141.59712
    ## 53       12   364       6   353  4.84e-61    213.0     56.83   138.88078
    ## 54       12   364       6   353  5.89e-61    213.0     56.83   138.68443
    ## 55       12   361      10   339  6.81e-61    214.0     56.25   138.53930
    ## 56       18   355       3   348  1.56e-60    212.0     55.56   137.71042
    ## 57       15   368      19   381  6.74e-60    211.0     55.64   136.24705
    ## 58       15   365       3   362  1.68e-59    209.0     55.29   135.33373
    ## 59       15   365       3   362  1.90e-59    209.0     55.29   135.21067
    ## 60       18   359     376   720  2.05e-59    218.0     58.40   135.13468
    ## 61       15   364       3   361  3.06e-59    208.0     55.44   134.73411
    ## 62       15   364       9   367  7.54e-58    205.0     55.17   131.52971
    ## 63       17   357       5   328  8.61e-58    203.0     54.86   131.39701
    ## 64       50   355      60   357  5.12e-56    199.0     56.82   127.31161
    ## 65       18   372      23   362  8.54e-56    199.0     52.78   126.80000
    ## 66       62   365      96   421  1.88e-55    200.0     54.52   126.01091
    ## 67       12   367      64   401  6.67e-52    189.0     53.30   117.83680
    ## 68       12   367      50   387  8.71e-52    188.0     53.30   117.56995
    ## 69       17   367      60   393  9.06e-52    188.0     53.20   117.53056
    ## 70       12   367      56   393  1.13e-51    188.0     53.30   117.30962
    ## 71       12   352      53   374  2.27e-51    186.0     53.74   116.61206
    ## 72       17   354      60   394  3.34e-51    186.0     51.66   116.22587
    ## 73       12   367      53   390  7.20e-51    186.0     53.02   115.45776
    ## 74       42   357      41   335  9.93e-51    183.0     55.00   115.13628
    ## 75       42   355      41   333  1.55e-50    183.0     55.03   114.69100
    ## 76       17   354      15   349  8.96e-50    181.0     47.86   112.93648
    ## 77       17   354       3   337  1.15e-49    180.0     47.86   112.68691
    ## 78       17   354       4   338  7.72e-49    178.0     47.59   110.78286
    ## 79       64   356      76   364  9.57e-49    178.0     56.31   110.56804
    ## 80       17   354       4   338  1.13e-48    177.0     47.59   110.40187
    ## 81       17   354       4   338  1.17e-48    177.0     47.59   110.36708
    ## 82       18   357      46   369  1.39e-48    178.0     53.16   110.19478
    ## 83       18   357      53   382  2.78e-47    175.0     53.03   107.19905
    ## 84       18   357       1   330  9.30e-47    171.0     51.87   105.99148
    ## 85       17   354       6   341  2.20e-46    171.0     52.10   105.13046
    ## 86       18   357      46   375  2.25e-46    177.0     52.84   105.10798
    ## 87       18   357      73   402  2.66e-46    172.0     51.87   104.94059
    ## 88       62   357     137   418  2.83e-46    173.0     54.36   104.87864
    ## 89       17   353      23   372  7.50e-46    170.0     50.14   103.90401
    ## 90       18   357      91   420  9.81e-46    171.0     53.47   103.63551
    ## 91       40   356     408   710  9.19e-45    174.0     54.80   101.39821
    ## 92       64   360     133   417  2.48e-44    167.0     54.61   100.40549
    ## 93       64   360      88   372  3.02e-44    166.0     54.61   100.20849
    ## 94       18   356       2   330  1.94e-40    154.0     52.14    91.44072
    ## 95      152   364     278   493  3.93e-36    145.0     60.71    81.52442
    ## 96       18   145      44   170  1.23e-07     56.6     47.83    15.91108
    ## 97       17   270       7   238  1.54e-32    127.0     52.14    73.25094
    ## 98       56   370      73   447  2.51e-31    129.0     46.97    70.45986
    ## 99      283   368       1    87  6.21e-21     89.4     67.82    46.52813
    ## 100     287   368       1    83  4.58e-20     87.8     68.67    44.53000
    ## 101     364   398       3    37  8.31e-07     48.5     77.14    14.00064
    ##     pdb.id    acc
    ## 1   3WPN_A 3WPN_A
    ## 2   3HQD_A 3HQD_A
    ## 3   1II6_A 1II6_A
    ## 4   4ZCA_A 4ZCA_A
    ## 5   4ZHI_A 4ZHI_A
    ## 6   4AP0_A 4AP0_A
    ## 7   1Q0B_A 1Q0B_A
    ## 8   4A1Z_A 4A1Z_A
    ## 9   4A28_A 4A28_A
    ## 10  1X88_A 1X88_A
    ## 11  4B7B_A 4B7B_A
    ## 12  5ZO7_A 5ZO7_A
    ## 13  4AQV_C 4AQV_C
    ## 14  3ZCW_A 3ZCW_A
    ## 15  2WBE_C 2WBE_C
    ## 16  5MM4_K 5MM4_K
    ## 17  5MM7_K 5MM7_K
    ## 18  5M5I_C 5M5I_C
    ## 19  5MLV_G 5MLV_G
    ## 20  6S8M_K 6S8M_K
    ## 21  3B6U_A 3B6U_A
    ## 22  3B6V_A 3B6V_A
    ## 23  1GOJ_A 1GOJ_A
    ## 24  2VVG_A 2VVG_A
    ## 25  3ZFC_A 3ZFC_A
    ## 26  3ZFD_A 3ZFD_A
    ## 27  6MLQ_C 6MLQ_C
    ## 28  4BN2_A 4BN2_A
    ## 29  2Y5W_A 2Y5W_A
    ## 30  6A20_A 6A20_A
    ## 31  6A1Z_A 6A1Z_A
    ## 32  4HNA_K 4HNA_K
    ## 33  3J8X_K 3J8X_K
    ## 34  4A14_A 4A14_A
    ## 35  1MKJ_A 1MKJ_A
    ## 36  4ATX_C 4ATX_C
    ## 37  4UXT_C 4UXT_C
    ## 38  3LRE_A 3LRE_A
    ## 39  3J6H_K 3J6H_K
    ## 40  5ZBS_B 5ZBS_B
    ## 41  4LNU_K 4LNU_K
    ## 42  5OAM_K 5OAM_K
    ## 43  5ZBR_B 5ZBR_B
    ## 44  1BG2_A 1BG2_A
    ## 45  2XT3_A 2XT3_A
    ## 46  5LT3_A 5LT3_A
    ## 47  5LT1_A 5LT1_A
    ## 48  3WRD_A 3WRD_A
    ## 49  5HNY_K 5HNY_K
    ## 50  1T5C_A 1T5C_A
    ## 51  5HNW_K 5HNW_K
    ## 52  1SDM_A 1SDM_A
    ## 53  5GSZ_A 5GSZ_A
    ## 54  5GSY_K 5GSY_K
    ## 55  3H4S_A 3H4S_A
    ## 56  3GBJ_A 3GBJ_A
    ## 57  1IA0_K 1IA0_K
    ## 58  1I6I_A 1I6I_A
    ## 59  1VFV_A 1VFV_A
    ## 60  4OZQ_A 4OZQ_A
    ## 61  1I5S_A 1I5S_A
    ## 62  4UXO_C 4UXO_C
    ## 63  5WDE_A 5WDE_A
    ## 64  3NWN_A 3NWN_A
    ## 65  6NJE_A 6NJE_A
    ## 66  2OWM_A 2OWM_A
    ## 67  2NCD_A 2NCD_A
    ## 68  1CZ7_A 1CZ7_A
    ## 69  3U06_A 3U06_A
    ## 70  5W3D_A 5W3D_A
    ## 71  3L1C_A 3L1C_A
    ## 72  4ETP_A 4ETP_A
    ## 73  1N6M_A 1N6M_A
    ## 74  3PXN_A 3PXN_A
    ## 75  3DC4_A 3DC4_A
    ## 76  1F9T_A 1F9T_A
    ## 77  3KAR_A 3KAR_A
    ## 78  1F9W_A 1F9W_A
    ## 79  4GKR_A 4GKR_A
    ## 80  1F9U_A 1F9U_A
    ## 81  1F9V_A 1F9V_A
    ## 82  4Y05_A 4Y05_A
    ## 83  2HEH_A 2HEH_A
    ## 84  3EDL_D 3EDL_D
    ## 85  3T0Q_A 3T0Q_A
    ## 86  5MIO_C 5MIO_C
    ## 87  1V8J_A 1V8J_A
    ## 88  5XJA_A 5XJA_A
    ## 89  5WDH_A 5WDH_A
    ## 90  2GRY_A 2GRY_A
    ## 91  4H1G_A 4H1G_A
    ## 92  6B0I_K 6B0I_K
    ## 93  3J2U_K 3J2U_K
    ## 94  1RY6_A 1RY6_A
    ## 95  5ND2_C 5ND2_C
    ## 96  5ND2_C 5ND2_C
    ## 97  2KIN_A 2KIN_A
    ## 98  5X3E_A 5X3E_A
    ## 99  2KIN_B 2KIN_B
    ## 100 3KIN_B 3KIN_B
    ## 101 5JV3_A 5JV3_A
    ## 
    ## $raw
    ##         queryid subjectids identity alignmentlength mismatches gapopens
    ## 1   Query_47107     3WPN_A   85.366             369         52        2
    ## 2   Query_47107     3HQD_A   85.366             369         52        2
    ## 3   Query_47107     1II6_A   85.326             368         52        2
    ## 4   Query_47107     4ZCA_A   85.095             369         53        2
    ## 5   Query_47107     4ZHI_A   85.095             369         53        2
    ## 6   Query_47107     4AP0_A   85.326             368         52        2
    ## 7   Query_47107     1Q0B_A   85.675             363         51        1
    ## 8   Query_47107     4A1Z_A   85.054             368         53        2
    ## 9   Query_47107     4A28_A   85.054             368         53        2
    ## 10  Query_47107     1X88_A   86.313             358         48        1
    ## 11  Query_47107     4B7B_A   84.783             368         54        2
    ## 12  Query_47107     5ZO7_A   86.686             353         46        1
    ## 13  Query_47107     4AQV_C   83.924             367         57        2
    ## 14  Query_47107     3ZCW_A   86.494             348         46        1
    ## 15  Query_47107     2WBE_C   56.383             376        152        6
    ## 16  Query_47107     5MM4_K   52.332             386        145        8
    ## 17  Query_47107     5MM7_K   52.030             394        149        9
    ## 18  Query_47107     5M5I_C   50.000             368        165        6
    ## 19  Query_47107     5MLV_G   50.000             368        165        6
    ## 20  Query_47107     6S8M_K   50.000             368        165        6
    ## 21  Query_47107     3B6U_A   46.832             363        175        8
    ## 22  Query_47107     3B6V_A   41.538             390        180       10
    ## 23  Query_47107     1GOJ_A   44.321             361        174        8
    ## 24  Query_47107     2VVG_A   42.297             357        180        8
    ## 25  Query_47107     3ZFC_A   43.094             362        182        9
    ## 26  Query_47107     3ZFD_A   43.094             362        182        9
    ## 27  Query_47107     6MLQ_C   41.274             361        182       10
    ## 28  Query_47107     4BN2_A   44.384             365        172       12
    ## 29  Query_47107     2Y5W_A   41.711             374        187       10
    ## 30  Query_47107     6A20_A   38.520             392        204       10
    ## 31  Query_47107     6A1Z_A   38.520             392        204       10
    ## 32  Query_47107     4HNA_K   42.254             355        176        8
    ## 33  Query_47107     3J8X_K   42.254             355        176        8
    ## 34  Query_47107     4A14_A   41.143             350        176       10
    ## 35  Query_47107     1MKJ_A   42.254             355        176        8
    ## 36  Query_47107     4ATX_C   42.254             355        176        8
    ## 37  Query_47107     4UXT_C   42.017             357        173       11
    ## 38  Query_47107     3LRE_A   39.669             363        181        9
    ## 39  Query_47107     3J6H_K   41.972             355        177        9
    ## 40  Query_47107     5ZBS_B   38.462             377        196        9
    ## 41  Query_47107     4LNU_K   42.151             344        170        8
    ## 42  Query_47107     5OAM_K   39.669             363        181        9
    ## 43  Query_47107     5ZBR_B   38.462             377        196        9
    ## 44  Query_47107     1BG2_A   42.151             344        170        8
    ## 45  Query_47107     2XT3_A   40.571             350        178       10
    ## 46  Query_47107     5LT3_A   41.860             344        171        8
    ## 47  Query_47107     5LT1_A   41.860             344        171        8
    ## 48  Query_47107     3WRD_A   41.880             351        175        9
    ## 49  Query_47107     5HNY_K   41.737             357        174       11
    ## 50  Query_47107     1T5C_A   39.943             353        185        9
    ## 51  Query_47107     5HNW_K   41.176             357        176       11
    ## 52  Query_47107     1SDM_A   37.726             387        205       12
    ## 53  Query_47107     5GSZ_A   37.705             366        197        9
    ## 54  Query_47107     5GSY_K   37.705             366        197        9
    ## 55  Query_47107     3H4S_A   38.920             352        191       10
    ## 56  Query_47107     3GBJ_A   38.333             360        186        9
    ## 57  Query_47107     1IA0_K   36.220             381        198       10
    ## 58  Query_47107     1I6I_A   35.979             378        197       10
    ## 59  Query_47107     1VFV_A   35.979             378        197       10
    ## 60  Query_47107     4OZQ_A   40.220             363        178       11
    ## 61  Query_47107     1I5S_A   36.074             377        196       10
    ## 62  Query_47107     4UXO_C   35.544             377        198       10
    ## 63  Query_47107     5WDE_A   38.857             350        179        9
    ## 64  Query_47107     3NWN_A   38.961             308        176        5
    ## 65  Query_47107     6NJE_A   35.556             360        207        8
    ## 66  Query_47107     2OWM_A   39.650             343        151       11
    ## 67  Query_47107     2NCD_A   35.440             364        201       12
    ## 68  Query_47107     1CZ7_A   35.440             364        201       12
    ## 69  Query_47107     3U06_A   35.376             359        199       11
    ## 70  Query_47107     5W3D_A   35.440             364        201       12
    ## 71  Query_47107     3L1C_A   35.632             348        191       11
    ## 72  Query_47107     4ETP_A   36.188             362        180       10
    ## 73  Query_47107     1N6M_A   35.165             364        202       12
    ## 74  Query_47107     3PXN_A   34.688             320        180        7
    ## 75  Query_47107     3DC4_A   34.906             318        178        7
    ## 76  Query_47107     1F9T_A   35.027             374        168       10
    ## 77  Query_47107     3KAR_A   35.027             374        168       10
    ## 78  Query_47107     1F9W_A   34.759             374        169       10
    ## 79  Query_47107     4GKR_A   36.570             309        160        8
    ## 80  Query_47107     1F9U_A   34.759             374        169       10
    ## 81  Query_47107     1F9V_A   34.759             374        169       10
    ## 82  Query_47107     4Y05_A   36.494             348        189       11
    ## 83  Query_47107     2HEH_A   35.447             347        200        9
    ## 84  Query_47107     3EDL_D   34.870             347        202        9
    ## 85  Query_47107     3T0Q_A   36.415             357        187       13
    ## 86  Query_47107     5MIO_C   35.795             352        192       11
    ## 87  Query_47107     1V8J_A   34.870             347        202        9
    ## 88  Query_47107     5XJA_A   36.242             298        172        6
    ## 89  Query_47107     5WDH_A   36.585             369        183       14
    ## 90  Query_47107     2GRY_A   34.393             346        205        7
    ## 91  Query_47107     4H1G_A   35.294             323        183        8
    ## 92  Query_47107     6B0I_K   36.513             304        167        9
    ## 93  Query_47107     3J2U_K   36.513             304        167        9
    ## 94  Query_47107     1RY6_A   33.903             351        198       11
    ## 95  Query_47107     5ND2_C   40.625             224        114        7
    ## 96  Query_47107     5ND2_C   26.087             138         81        4
    ## 97  Query_47107     2KIN_A   37.743             257        132        8
    ## 98  Query_47107     5X3E_A   29.551             379        199       11
    ## 99  Query_47107     2KIN_B   55.172              87         38        1
    ## 100 Query_47107     3KIN_B   55.422              83         36        1
    ## 101 Query_47107     5JV3_A   71.429              35         10        0
    ##     q.start q.end s.start s.end    evalue bitscore positives
    ## 1         1   367       7   375  0.00e+00    677.0     94.58
    ## 2         1   367       1   369  0.00e+00    676.0     94.58
    ## 3         1   366       1   368  0.00e+00    675.0     94.57
    ## 4         1   367       1   369  0.00e+00    675.0     94.58
    ## 5         1   367       1   369  0.00e+00    674.0     94.31
    ## 6         1   366       3   370  0.00e+00    674.0     94.57
    ## 7         5   366       5   367  0.00e+00    674.0     95.04
    ## 8         1   366       1   368  0.00e+00    672.0     94.29
    ## 9         1   366       1   368  0.00e+00    672.0     94.29
    ## 10       10   366       2   359  0.00e+00    671.0     95.81
    ## 11        1   366       1   368  0.00e+00    668.0     94.02
    ## 12       16   367       8   360  0.00e+00    664.0     95.75
    ## 13        1   365       1   367  0.00e+00    658.0     93.19
    ## 14       15   361       1   348  0.00e+00    654.0     95.69
    ## 15        1   369       3   373 3.21e-136    416.0     74.73
    ## 16       17   368       5   385 4.11e-113    356.0     67.88
    ## 17       10   368      69   457 4.61e-113    358.0     67.51
    ## 18       12   365       4   366 1.25e-106    338.0     67.93
    ## 19       12   365      67   429 2.39e-105    337.0     67.93
    ## 20       12   365      73   435 6.10e-105    336.0     67.93
    ## 21       18   376      23   371  4.11e-95    307.0     63.91
    ## 22       18   376      23   395  1.13e-84    280.0     59.23
    ## 23       17   372       7   345  7.03e-82    271.0     61.22
    ## 24       17   365       5   343  3.60e-81    269.0     61.06
    ## 25        9   364       1   344  2.67e-79    264.0     61.33
    ## 26        9   364       1   344  2.76e-79    264.0     61.33
    ## 27       18   366      17   359  1.35e-69    239.0     59.56
    ## 28       15   365      10   357  5.33e-69    236.0     56.71
    ## 29       18   385      13   361  3.40e-68    234.0     57.49
    ## 30       18   387       6   382  4.49e-68    236.0     56.38
    ## 31       18   387       6   382  1.16e-67    233.0     56.38
    ## 32       17   368       8   336  3.99e-67    230.0     55.21
    ## 33       17   368       8   336  8.63e-67    229.0     55.21
    ## 34       18   355      13   344  1.52e-66    228.0     59.43
    ## 35       17   368       8   336  2.64e-66    228.0     55.21
    ## 36       17   368       8   336  4.99e-66    227.0     54.93
    ## 37       18   368      10   338  1.67e-65    225.0     56.30
    ## 38       16   357      10   355  5.79e-65    224.0     56.75
    ## 39       18   368       9   338  6.04e-65    224.0     55.77
    ## 40       18   372       4   366  6.36e-65    225.0     55.70
    ## 41       17   357       8   325  8.37e-65    223.0     55.23
    ## 42       16   357      13   358  1.10e-64    224.0     56.75
    ## 43       18   372       4   366  1.29e-64    224.0     55.70
    ## 44       17   357       8   325  2.93e-64    221.0     55.23
    ## 45       18   355      13   344  3.34e-64    222.0     59.14
    ## 46       17   357       8   325  4.57e-64    221.0     54.94
    ## 47       17   357       8   325  5.25e-64    220.0     54.94
    ## 48       18   364       9   334  7.15e-64    221.0     55.56
    ## 49       17   367      24   352  1.10e-63    221.0     56.30
    ## 50       18   365       6   336  4.56e-63    219.0     55.52
    ## 51       17   367      24   352  2.07e-62    218.0     56.02
    ## 52       12   395       2   355  3.20e-62    217.0     56.59
    ## 53       12   364       6   353  4.84e-61    213.0     56.83
    ## 54       12   364       6   353  5.89e-61    213.0     56.83
    ## 55       12   361      10   339  6.81e-61    214.0     56.25
    ## 56       18   355       3   348  1.56e-60    212.0     55.56
    ## 57       15   368      19   381  6.74e-60    211.0     55.64
    ## 58       15   365       3   362  1.68e-59    209.0     55.29
    ## 59       15   365       3   362  1.90e-59    209.0     55.29
    ## 60       18   359     376   720  2.05e-59    218.0     58.40
    ## 61       15   364       3   361  3.06e-59    208.0     55.44
    ## 62       15   364       9   367  7.54e-58    205.0     55.17
    ## 63       17   357       5   328  8.61e-58    203.0     54.86
    ## 64       50   355      60   357  5.12e-56    199.0     56.82
    ## 65       18   372      23   362  8.54e-56    199.0     52.78
    ## 66       62   365      96   421  1.88e-55    200.0     54.52
    ## 67       12   367      64   401  6.67e-52    189.0     53.30
    ## 68       12   367      50   387  8.71e-52    188.0     53.30
    ## 69       17   367      60   393  9.06e-52    188.0     53.20
    ## 70       12   367      56   393  1.13e-51    188.0     53.30
    ## 71       12   352      53   374  2.27e-51    186.0     53.74
    ## 72       17   354      60   394  3.34e-51    186.0     51.66
    ## 73       12   367      53   390  7.20e-51    186.0     53.02
    ## 74       42   357      41   335  9.93e-51    183.0     55.00
    ## 75       42   355      41   333  1.55e-50    183.0     55.03
    ## 76       17   354      15   349  8.96e-50    181.0     47.86
    ## 77       17   354       3   337  1.15e-49    180.0     47.86
    ## 78       17   354       4   338  7.72e-49    178.0     47.59
    ## 79       64   356      76   364  9.57e-49    178.0     56.31
    ## 80       17   354       4   338  1.13e-48    177.0     47.59
    ## 81       17   354       4   338  1.17e-48    177.0     47.59
    ## 82       18   357      46   369  1.39e-48    178.0     53.16
    ## 83       18   357      53   382  2.78e-47    175.0     53.03
    ## 84       18   357       1   330  9.30e-47    171.0     51.87
    ## 85       17   354       6   341  2.20e-46    171.0     52.10
    ## 86       18   357      46   375  2.25e-46    177.0     52.84
    ## 87       18   357      73   402  2.66e-46    172.0     51.87
    ## 88       62   357     137   418  2.83e-46    173.0     54.36
    ## 89       17   353      23   372  7.50e-46    170.0     50.14
    ## 90       18   357      91   420  9.81e-46    171.0     53.47
    ## 91       40   356     408   710  9.19e-45    174.0     54.80
    ## 92       64   360     133   417  2.48e-44    167.0     54.61
    ## 93       64   360      88   372  3.02e-44    166.0     54.61
    ## 94       18   356       2   330  1.94e-40    154.0     52.14
    ## 95      152   364     278   493  3.93e-36    145.0     60.71
    ## 96       18   145      44   170  1.23e-07     56.6     47.83
    ## 97       17   270       7   238  1.54e-32    127.0     52.14
    ## 98       56   370      73   447  2.51e-31    129.0     46.97
    ## 99      283   368       1    87  6.21e-21     89.4     67.82
    ## 100     287   368       1    83  4.58e-20     87.8     68.67
    ## 101     364   398       3    37  8.31e-07     48.5     77.14
    ## 
    ## $url
    ##                                                                                                                                                        YWWTRSBY014 
    ## "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&RESULTS_FILE=on&FORMAT_TYPE=CSV&ALIGNMENTS=20000&RID=YWWTRSBY014" 
    ## 
    ## attr(,"class")
    ## [1] "blast"

``` r
#write.csv(results,"blastresults.txt")
```
