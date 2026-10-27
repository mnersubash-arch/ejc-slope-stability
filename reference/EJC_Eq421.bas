Attribute VB_Name = "EJC_Eq421"
'=====================================================================================================
' EJC_Eq421 (GENG5512 package v14; library revision 6, corrected scan boundary; v14a validation-sheet text fix)
' Instantaneous friction angle phi'i from eq. 4.21 (Fu & Liao 2010 / Kumar 1998 C1 with the rock-mass s)
' for a principal-stress closure  sigma2' = pC*sigma1' + qC*sigma3', and the acceptance checks of
' methodology Section 4.2. Strength library called by the primary VBA assembly and the worksheet QA.
'   Closure 1 (thesis):  pC = -1, qC = 1/nu     Closure 2 (Moscarda): pC = qC = nu     Closure 3 (2D HB): pC = 0, qC = 1
' Units: sigma_n' and c'i in kPa; sigma_ci in MPa; angles in degrees (inputs/outputs), radians internally.
' Solver: finite scan from max(phi_b + 0.001 deg, 0.01 deg) to 89.9999 deg (1800 intervals), sign changes between
'         DEFINED points with a preliminary pole filter, then Newton-Raphson with the ANALYTICAL derivative of
'         g(phi) = sin(phi) - (M-1)/(M+1), kept inside the bracket. Converged only when the UNCLIPPED angular
'         residual |phi - asin((M-1)/(M+1))| <= tol (radians) AND the last step <= tol. Otherwise an error value.
' Adopted branch (same as the workbook solver and admissibility table): 0 <= phi'i < 90 deg and tangent angle
' asin((M-1)/(M+1)) >= -tol (i.e. M >= 1 up to roundoff). States outside it are REJECTED: DOMAIN / BRANCH.
' EJC_NBRACKETS counts brackets OBSERVED by the finite scan; it is not proof that every root has been found.
' Error values: #VALUE! invalid input (incl. rootIndex < 1, tol <= 0), #N/A no (further) root bracket, #NUM! domain failure / not converged.
' Algorithm mirrored in Python (EJC_Eq421_port.py, 22 shared fixtures; assembly/fast-HB paths have separate tests). In Excel: Debug > Compile
' VBAProject, then run EJC_RunValidation and keep the 'VBA Validation' sheet as the test record.
'=====================================================================================================
Option Explicit
Private Const PI_ As Double = 3.14159265358979
Private Const HI_DEG As Double = 89.9999
Private Const NSCAN As Long = 1800

Private Function ValidInputs(ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, ByVal tolA As Double) As Boolean
    ValidInputs = (sciV > 0 And mbV > 0 And sV >= 0 And aV > 0 And aV <= 1 And tolA > 0 And tolA <= 0.000001)
End Function

Private Function GState(ByVal phi As Double, ByVal sn As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, _
                        ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, ByRef gOut As Double, ByRef dgOut As Double, _
                        ByRef mOut As Double) As Boolean
    ' g(phi) = sin(phi) - 1 + 2D/(N+D) (= sin(phi) - (M-1)/(M+1), M = N/D) and its exact derivative
    Dim cP As Double, cR As Double, cS As Double, sinP As Double, cosP As Double, bTerm As Double
    Dim hT As Double, dhT As Double, sig1 As Double, sig3 As Double, dsig1 As Double, dsig3 As Double, hbQ As Double, dhbQ As Double
    Dim nT As Double, dnT As Double, dT As Double, ddT As Double
    GState = False
    cP = pC * pC - pC + 1: cR = 2 * pC * qC - pC - qC - 1: cS = qC * qC - qC + 1
    sinP = Sin(phi): cosP = Cos(phi): bTerm = mbV * sn / sciV + sV
    If bTerm <= 0 Then Exit Function
    hT = (bTerm / (1 + sinP / aV)) ^ aV: dhT = -hT * cosP / (1 + sinP / aV)
    sig3 = sn - sciV / 2 * hT * (1 - sinP): dsig3 = -sciV / 2 * (dhT * (1 - sinP) - hT * cosP)
    sig1 = sn + sciV / 2 * hT * (1 + sinP): dsig1 = sciV / 2 * (dhT * (1 + sinP) + hT * cosP)
    hbQ = mbV * sig3 / sciV + sV: dhbQ = mbV * dsig3 / sciV
    If hbQ <= 0 Then Exit Function                                   ' beyond the Hoek-Brown tensile limit
    nT = 2 * aV * mbV * sciV * hbQ ^ (2 * aV - 1) - cR * sig1 - 2 * cS * sig3
    dnT = 2 * aV * mbV * sciV * (2 * aV - 1) * hbQ ^ (2 * aV - 2) * dhbQ - cR * dsig1 - 2 * cS * dsig3
    dT = 2 * cP * sig1 + cR * sig3: ddT = 2 * cP * dsig1 + cR * dsig3
    If nT + dT = 0 Or dT = 0 Then Exit Function
    gOut = sinP - 1 + 2 * dT / (nT + dT)
    dgOut = cosP - 2 * (dnT * dT - nT * ddT) / (nT + dT) ^ 2
    mOut = nT / dT
    GState = True
End Function

Private Function AngRes(ByVal phi As Double, ByVal sn As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, _
                        ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, ByRef okOut As Boolean) As Double
    ' UNCLIPPED eq. 4.19 residual phi - asin((M-1)/(M+1)); okOut = False if M+1 = 0 or |argument| > 1
    Dim gV As Double, dgV As Double, mV As Double, zV As Double
    okOut = False
    If Not GState(phi, sn, sciV, mbV, sV, aV, pC, qC, gV, dgV, mV) Then Exit Function
    If mV + 1 = 0 Then Exit Function
    zV = (mV - 1) / (mV + 1)
    If Abs(zV) > 1 Then Exit Function
    AngRes = phi - Application.WorksheetFunction.Asin(zV)
    okOut = True
End Function

Private Function PhiB(ByVal sn As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, _
                      ByRef okOut As Boolean) As Double
    ' angle above which the eq. 4.21 state is inside the tensile limit (Q > 0): closed form (a = 0.5) + Newton
    Dim tLim As Double, kT As Double, zT As Double, pV As Double, fV As Double, fpV As Double, it As Integer
    okOut = False: tLim = sV * sciV / mbV: kT = sn + tLim
    If kT <= 0 Then Exit Function
    zT = 2 * kT / (sciV * (mbV * kT / sciV) ^ aV)
    okOut = True
    If zT >= 1 Then PhiB = 0: Exit Function
    pV = (1 + zT * zT) - Sqr((1 + zT * zT) ^ 2 - (1 - zT * zT))
    For it = 1 To 8
        fV = (1 - pV) * (1 + pV / aV) ^ (-aV) - zT
        fpV = -(1 + pV / aV) ^ (-aV) * (1 + (1 - pV) / (1 + pV / aV))
        If fpV = 0 Then Exit For
        pV = pV - fV / fpV
        If pV < 0 Then pV = 0
        If pV > 1 - 0.000000000000001 Then pV = 1 - 0.000000000000001
    Next it
    PhiB = Application.WorksheetFunction.Asin(pV)
End Function

Private Function FindBracket(ByVal sn As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, _
                             ByVal pC As Double, ByVal qC As Double, ByVal wantIndex As Long, ByRef loOut As Double, ByRef hiOut As Double) As Long
    ' returns observed candidate brackets after a preliminary pole filter; final angular/state checks still decide admissibility.
    Dim pb As Double, okB As Boolean, x0 As Double, x1 As Double, k As Long, x As Double, gV As Double, dgV As Double, mV As Double
    Dim havePrev As Boolean, xPrev As Double, gPrev As Double, cnt As Long
    Dim lo As Double, hi As Double, glo As Double, ghi As Double, xm As Double, gm As Double, dgm As Double, mm As Double, j As Integer, isPole As Boolean
    cnt = 0: FindBracket = 0
    pb = PhiB(sn, sciV, mbV, sV, aV, okB)
    If Not okB Then Exit Function
    ' A point only 1E-12 rad above the tensile boundary can create a spurious first bracket.
    ' Keep the candidate-selection policy; move the finite scan start 0.001 degree inside the domain.
    x0 = pb + 0.001 * PI_ / 180
    If x0 < 0.01 * PI_ / 180 Then x0 = 0.01 * PI_ / 180
    x1 = HI_DEG * PI_ / 180
    If x0 >= x1 Then Exit Function
    havePrev = False
    For k = 0 To NSCAN
        x = x0 + (x1 - x0) * k / NSCAN
        If GState(x, sn, sciV, mbV, sV, aV, pC, qC, gV, dgV, mV) Then
            If havePrev Then
                If gPrev * gV <= 0 Then
                    lo = xPrev: hi = x: glo = gPrev: ghi = gV: isPole = False
                    For j = 1 To 30                                  ' bisect; |g| that stays large means a pole
                        xm = (lo + hi) / 2
                        If Not GState(xm, sn, sciV, mbV, sV, aV, pC, qC, gm, dgm, mm) Then isPole = True: Exit For
                        If glo * gm <= 0 Then
                            hi = xm: ghi = gm
                        Else
                            lo = xm: glo = gm
                        End If
                    Next j
                    If Not isPole Then
                        If Abs(glo) > 0.001 And Abs(ghi) > 0.001 Then isPole = True
                    End If
                    If Not isPole Then
                        cnt = cnt + 1
                        If cnt = wantIndex Then loOut = xPrev: hiOut = x
                    End If
                End If
            End If
            havePrev = True: xPrev = x: gPrev = gV
        Else
            havePrev = False
        End If
    Next k
    FindBracket = cnt
End Function

Public Function EJC_PHI(ByVal sn_kPa As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, _
                        ByVal pC As Double, ByVal qC As Double, Optional ByVal tolA As Double = 0.0000000001, _
                        Optional ByVal rootIndex As Long = 1) As Variant
    ' CANDIDATE root only (call EJC_ACCEPT / EJC_CHECK for acceptance). Degrees, or an error value.
    Dim sn As Double, lo As Double, hi As Double, nB As Long, x As Double, xn As Double, gV As Double, dgV As Double, mV As Double
    Dim glo As Double, dgl As Double, ml As Double, it As Integer, stepV As Double, rV As Double, okR As Boolean
    If Not ValidInputs(sciV, mbV, sV, aV, tolA) Or rootIndex < 1 Then EJC_PHI = CVErr(xlErrValue): Exit Function
    sn = sn_kPa / 1000#
    nB = FindBracket(sn, sciV, mbV, sV, aV, pC, qC, rootIndex, lo, hi)
    If nB < rootIndex Then EJC_PHI = CVErr(xlErrNA): Exit Function
    If Not GState(lo, sn, sciV, mbV, sV, aV, pC, qC, glo, dgl, ml) Then EJC_PHI = CVErr(xlErrNum): Exit Function
    x = (lo + hi) / 2
    For it = 1 To 200
        If Not GState(x, sn, sciV, mbV, sV, aV, pC, qC, gV, dgV, mV) Then EJC_PHI = CVErr(xlErrNum): Exit Function
        If glo * gV <= 0 Then
            hi = x
        Else
            lo = x: glo = gV
        End If
        If dgV <> 0 Then xn = x - gV / dgV Else xn = (lo + hi) / 2
        If xn <= lo Or xn >= hi Then xn = (lo + hi) / 2
        stepV = Abs(xn - x): x = xn
        rV = AngRes(x, sn, sciV, mbV, sV, aV, pC, qC, okR)
        If okR Then
            If Abs(rV) <= tolA And stepV <= tolA Then EJC_PHI = x * 180 / PI_: Exit Function
        End If
    Next it
    EJC_PHI = CVErr(xlErrNum)                                         ' iterations exhausted: not converged
End Function

Public Function EJC_FASTPHI(ByVal sn_kPa As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, _
                           ByVal pC As Double, ByVal qC As Double, ByVal guessDeg As Double, Optional ByVal tolA As Double = 0.000000000001) As Variant
    ' Same eq. 4.21 and analytical derivative. The conventional-HB branch has one tangent root.
    ' Its tensile-domain bound provides a single bracket; a warm start changes only the numerical path.
    ' Other closures retain the scanned FIRST candidate root, with no root switching or fallback strength.
    Dim sn As Double, lo As Double, hi As Double, x As Double, xn As Double, pb As Double
    Dim gv As Double, dgv As Double, mv As Double, glo As Double, ghi As Double, angleResidual As Double
    Dim okB As Boolean, okA As Boolean, it As Long, stepV As Double
    If Not ValidInputs(sciV, mbV, sV, aV, tolA) Then EJC_FASTPHI = CVErr(xlErrValue): Exit Function
    If pC <> 0 Or qC <> 1 Then
        EJC_FASTPHI = EJC_PHI(sn_kPa, sciV, mbV, sV, aV, pC, qC, tolA): Exit Function
    End If
    sn = sn_kPa / 1000
    pb = PhiB(sn, sciV, mbV, sV, aV, okB)
    If Not okB Then EJC_FASTPHI = CVErr(xlErrNum): Exit Function
    lo = pb + 0.000000001: hi = HI_DEG * PI_ / 180
    If lo >= hi Then EJC_FASTPHI = CVErr(xlErrNum): Exit Function
    If Not GState(lo, sn, sciV, mbV, sV, aV, pC, qC, glo, dgv, mv) Then GoTo FullScan
    If Not GState(hi, sn, sciV, mbV, sV, aV, pC, qC, ghi, dgv, mv) Then GoTo FullScan
    If glo * ghi > 0 Then GoTo FullScan
    x = (lo + hi) / 2
    If guessDeg * PI_ / 180 > lo And guessDeg * PI_ / 180 < hi Then x = guessDeg * PI_ / 180
    For it = 1 To 200
        If Not GState(x, sn, sciV, mbV, sV, aV, pC, qC, gv, dgv, mv) Then GoTo FullScan
        angleResidual = AngRes(x, sn, sciV, mbV, sV, aV, pC, qC, okA)
        If okA Then
            If dgv <> 0 Then
                If Abs(angleResidual) <= tolA And Abs(gv / dgv) <= tolA Then EJC_FASTPHI = x * 180 / PI_: Exit Function
            End If
        End If
        If glo * gv <= 0 Then hi = x Else lo = x: glo = gv
        If dgv <> 0 Then xn = x - gv / dgv Else xn = (lo + hi) / 2
        If xn <= lo Or xn >= hi Then xn = (lo + hi) / 2
        stepV = Abs(xn - x): x = xn
        angleResidual = AngRes(x, sn, sciV, mbV, sV, aV, pC, qC, okA)
        If okA Then
            If Abs(angleResidual) <= tolA And stepV <= tolA Then EJC_FASTPHI = x * 180 / PI_: Exit Function
        End If
    Next it
FullScan:
    EJC_FASTPHI = EJC_PHI(sn_kPa, sciV, mbV, sV, aV, pC, qC, tolA)
End Function

Public Function EJC_NBRACKETS(ByVal sn_kPa As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, _
                              ByVal pC As Double, ByVal qC As Double) As Variant
    Dim lo As Double, hi As Double
    If Not ValidInputs(sciV, mbV, sV, aV, 0.0000000001) Then EJC_NBRACKETS = CVErr(xlErrValue): Exit Function
    EJC_NBRACKETS = FindBracket(sn_kPa / 1000#, sciV, mbV, sV, aV, pC, qC, 0, lo, hi)
End Function

Private Function DomainFail() As Variant
    Dim res(1 To 8) As Variant, i As Integer
    For i = 1 To 8: res(i) = CVErr(xlErrNum): Next i
    res(4) = 0: res(7) = 0: res(8) = 0
    DomainFail = res
End Function

Private Function StateChecks(ByVal phi As Double, ByVal sig1 As Double, ByVal sig2 As Double, ByVal sig3 As Double, ByVal sciV As Double, _
                             ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, _
                             ByVal tolA As Double) As Variant
    ' 1 x 8: [angular residual (rad), original unsquared EJC residual (MPa), closure residual (MPa), ordering (1 = OK),
    '         sigma1' (MPa), sigma3' (MPa), domain (1 = derivative/radicand conditions valid), branch (1 = adopted branch)] - no clipping
    Dim res(1 To 8) As Variant, cP As Double, cR As Double, cS As Double, hbQ As Double, denM As Double, mV As Double, zV As Double, radV As Double
    Dim tAng As Double
    cP = pC * pC - pC + 1: cR = 2 * pC * qC - pC - qC - 1: cS = qC * qC - qC + 1
    res(1) = CVErr(xlErrNum): res(2) = CVErr(xlErrNum): res(3) = sig2 - (pC * sig1 + qC * sig3)
    res(4) = IIf(sig1 >= sig2 - 0.000000001 And sig2 >= sig3 - 0.000000001, 1, 0)
    res(5) = sig1: res(6) = sig3: res(7) = 0: res(8) = 0
    hbQ = mbV * sig3 / sciV + sV
    If hbQ > 0 Then
        denM = 2 * cP * sig1 + cR * sig3
        If denM <> 0 Then
            mV = (2 * aV * mbV * sciV * hbQ ^ (2 * aV - 1) - cR * sig1 - 2 * cS * sig3) / denM
            If mV + 1 <> 0 Then
                zV = (mV - 1) / (mV + 1)
                radV = 4 * sciV ^ 2 * hbQ ^ (2 * aV) - 3 * (sig2 - sig3) ^ 2
                If Abs(zV) <= 1 And radV >= 0 Then
                    tAng = Application.WorksheetFunction.Asin(zV)
                    res(1) = phi - tAng
                    res(2) = sig1 - (sig2 + sig3 + Sqr(radV)) / 2
                    res(7) = 1
                    If phi >= 0 And phi < PI_ / 2 And tAng >= -tolA Then res(8) = 1
                End If
            End If
        End If
    End If
    StateChecks = res
End Function

Public Function EJC_CHECK(ByVal sn_kPa As Double, ByVal phi_deg As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, _
                          ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, Optional ByVal tolA As Double = 0.0000000001) As Variant
    ' checks of the eq. 4.21 (C1) CANDIDATE state at (sigma_n', phi'i). Enter as an array formula over 8 cells.
    Dim sn As Double, phi As Double, hT As Double, sig1 As Double, sig3 As Double, bTerm As Double, dTerm As Double
    If Not ValidInputs(sciV, mbV, sV, aV, tolA) Then EJC_CHECK = DomainFail(): Exit Function
    If phi_deg < 0 Or phi_deg >= 90 Then EJC_CHECK = DomainFail(): Exit Function          ' adopted branch 0 <= phi < 90 deg
    sn = sn_kPa / 1000#: phi = phi_deg * PI_ / 180
    bTerm = mbV * sn / sciV + sV: dTerm = 1 + Sin(phi) / aV
    If bTerm <= 0 Or dTerm <= 0 Then EJC_CHECK = DomainFail(): Exit Function             ' guard before the fractional power
    hT = (bTerm / dTerm) ^ aV
    sig3 = sn - sciV / 2 * hT * (1 - Sin(phi)): sig1 = sn + sciV / 2 * hT * (1 + Sin(phi))
    EJC_CHECK = StateChecks(phi, sig1, pC * sig1 + qC * sig3, sig3, sciV, mbV, sV, aV, pC, qC, tolA)
End Function

Public Function EJC_CHECK_STATE(ByVal sn_kPa As Double, ByVal c_kPa As Double, ByVal phi_deg As Double, ByVal sciV As Double, _
                                ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, _
                                Optional ByVal tolA As Double = 0.0000000001) As Variant
    ' checks of the FORCE-USED state (sigma_n', c'i, phi'i): tau = c + sigma_n tan(phi), centre sigma_n + tau tan(phi), radius tau/cos(phi)
    Dim phi As Double, tauV As Double, cen As Double, radius As Double, sig1 As Double, sig3 As Double
    If Not ValidInputs(sciV, mbV, sV, aV, tolA) Then EJC_CHECK_STATE = DomainFail(): Exit Function
    If phi_deg < 0 Or phi_deg >= 90 Then EJC_CHECK_STATE = DomainFail(): Exit Function
    phi = phi_deg * PI_ / 180
    tauV = c_kPa + sn_kPa * Tan(phi): cen = sn_kPa + tauV * Tan(phi): radius = tauV / Cos(phi)
    sig1 = (cen + radius) / 1000#: sig3 = (cen - radius) / 1000#
    EJC_CHECK_STATE = StateChecks(phi, sig1, pC * sig1 + qC * sig3, sig3, sciV, mbV, sV, aV, pC, qC, tolA)
End Function

Private Function AcceptText(ByVal chk As Variant, ByVal tolA As Double, ByVal tolE As Double) As String
    If Not (tolA > 0 And tolA <= 0.000001 And tolE > 0 And tolE <= 0.0001) Then AcceptText = "REJECTED: INPUT": Exit Function
    If Not IsArray(chk) Then AcceptText = "REJECTED: INPUT/DOMAIN": Exit Function
    If chk(7) <> 1 Then AcceptText = "REJECTED: DOMAIN": Exit Function
    If chk(8) <> 1 Then AcceptText = "REJECTED: BRANCH": Exit Function
    If Abs(chk(1)) > tolA Then AcceptText = "REJECTED: TANGENCY": Exit Function
    If Abs(chk(2)) > tolE Then AcceptText = "REJECTED: EJC": Exit Function
    If Abs(chk(3)) > tolE Then AcceptText = "REJECTED: CLOSURE": Exit Function
    If chk(4) <> 1 Then AcceptText = "REJECTED: ORDERING": Exit Function
    AcceptText = "ADMISSIBLE"
End Function

Public Function EJC_ACCEPT(ByVal sn_kPa As Double, ByVal phi_deg As Double, ByVal sciV As Double, ByVal mbV As Double, ByVal sV As Double, _
                           ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, Optional ByVal tolA As Double = 0.0000000001, _
                           Optional ByVal tolE As Double = 0.00000001) As String
    EJC_ACCEPT = AcceptText(EJC_CHECK(sn_kPa, phi_deg, sciV, mbV, sV, aV, pC, qC, tolA), tolA, tolE)
End Function

Public Function EJC_ACCEPT_STATE(ByVal sn_kPa As Double, ByVal c_kPa As Double, ByVal phi_deg As Double, ByVal sciV As Double, _
                                 ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, _
                                 Optional ByVal tolA As Double = 0.0000000001, Optional ByVal tolE As Double = 0.00000001) As String
    EJC_ACCEPT_STATE = AcceptText(EJC_CHECK_STATE(sn_kPa, c_kPa, phi_deg, sciV, mbV, sV, aV, pC, qC, tolA), tolA, tolE)
End Function

Private Function ErrText(ByVal v As Variant) As String
    ' readable text for a result that may be an Excel error value (comparison per Microsoft 'Cell Error Values')
    If Not IsError(v) Then ErrText = CStr(v): Exit Function
    Select Case v
        Case CVErr(xlErrNA): ErrText = "#N/A"
        Case CVErr(xlErrNum): ErrText = "#NUM!"
        Case CVErr(xlErrValue): ErrText = "#VALUE!"
        Case Else: ErrText = "error"
    End Select
End Function

Private Function IsErrCode(ByVal v As Variant, ByVal code As Long) As Boolean
    IsErrCode = False
    ' CVErr has subtype Error; use documented explicit conversion before comparing codes.
    If IsError(v) Then IsErrCode = (CLng(v) = code)
End Function

Private Function NumTxt(ByVal v As Variant) As String
    If IsError(v) Then NumTxt = ErrText(v) Else NumTxt = Format(v, "0.000E+00")
End Function

Private Sub RecordTest(ByRef ws As Worksheet, ByVal rowN As Long, ByVal nm As String, ByVal ex As String, ByVal act As String, _
                       ByVal chk As Variant, ByVal ok As Boolean)
    ' v14a: the Expected label is written as text. Without the apostrophe Excel parses the label "#VALUE!"
    ' (T6, T10, T11, T13, T14) into a real #VALUE! error constant, which the workbook QA then counts.
    ws.Cells(rowN, 1).Value = "'" & nm: ws.Cells(rowN, 2).Value = "'" & ex: ws.Cells(rowN, 3).Value = "'" & act
    If IsArray(chk) Then
        ws.Cells(rowN, 4).Value = "'" & NumTxt(chk(1)): ws.Cells(rowN, 5).Value = "'" & NumTxt(chk(2))
        ws.Cells(rowN, 6).Value = AcceptText(chk, 0.0000000001, 0.00000001)
    Else
        ws.Cells(rowN, 4).Value = "-": ws.Cells(rowN, 5).Value = "-": ws.Cells(rowN, 6).Value = "-"
    End If
    ws.Cells(rowN, 7).Value = IIf(ok, "PASS", "FAIL")
End Sub

Public Sub EJC_RunValidation()
    ' Writes 22 shared fixtures. G26 remains the stable summary gate; T22 is recorded at row 28.
    Dim ws As Worksheet, mbC As Double, sC As Double, ph As Variant, ok As Boolean, chk As Variant, lbl As String
    Dim cV As Double, phT1 As Double, nPass As Long, i As Long, hit As Boolean, rw As Long, emptyChk As Variant
    Const T7SN As Double = 13265.01174607056
    mbC = 10 * Exp((42 - 100) / 28): sC = Exp((42 - 100) / 9)
    On Error Resume Next
    Set ws = ThisWorkbook.Worksheets("VBA Validation")
    On Error GoTo 0
    If ws Is Nothing Then
        Set ws = ThisWorkbook.Worksheets.Add(After:=ThisWorkbook.Worksheets(ThisWorkbook.Worksheets.Count))
        ws.Name = "VBA Validation"
    End If
    ws.Cells.Clear
    ws.Range("A1").Value = "v14 / EJC_Eq421 revision 6 validation - Excel " & Application.Version & " - " & Format(Now, "yyyy-mm-dd hh:nn")
    ws.Range("A3").Value = "Test": ws.Range("B3").Value = "Expected": ws.Range("C3").Value = "Returned"
    ws.Range("D3").Value = "Angular residual (rad)": ws.Range("E3").Value = "EJC residual (MPa)": ws.Range("F3").Value = "Acceptance label": ws.Range("G3").Value = "Result"
    emptyChk = Empty: rw = 4: nPass = 0
    ' T1
    ph = EJC_PHI(400, 77.7, mbC, sC, 0.5, 0, 1): ok = Not IsError(ph): chk = emptyChk
    If ok Then chk = EJC_CHECK(400, ph, 77.7, mbC, sC, 0.5, 0, 1): phT1 = ph: ok = (Abs(ph - 58.448348604037) < 0.000000001 And AcceptText(chk, 0.0000000001, 0.00000001) = "ADMISSIBLE")
    RecordTest ws, rw, "T1 closure 3, Chile, 400 kPa", "58.448348604037, ADMISSIBLE", ErrText(ph), chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T2
    ph = EJC_PHI(5740.3, 120, 3.5, 0.002, 0.505, -1, 20): ok = Not IsError(ph): chk = emptyChk
    If ok Then chk = EJC_CHECK(5740.3, ph, 120, 3.5, 0.002, 0.505, -1, 20): ok = (Abs(ph - 53.3178012) < 0.000001 And AcceptText(chk, 0.0000000001, 0.00000001) = "REJECTED: EJC")
    RecordTest ws, rw, "T2 closure 1, verification point, first root", "53.3178012, REJECTED: EJC", ErrText(ph), chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T3 (honest outcome: a converged admissible root OR a controlled error)
    ph = EJC_PHI(-98.00411947453206, 77.7, mbC, sC, 0.5, 0, 1): chk = emptyChk
    If IsError(ph) Then
        ok = True
    Else
        chk = EJC_CHECK(-98.00411947453206, ph, 77.7, mbC, sC, 0.5, 0, 1): ok = (AcceptText(chk, 0.0000000001, 0.00000001) = "ADMISSIBLE")
    End If
    RecordTest ws, rw, "T3 near tensile limit (-98.0041 kPa)", "about 89.6614 admissible, or a controlled error", ErrText(ph), chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T4 (asserted near-limit success)
    ph = EJC_PHI(-98.0041250683277, 77.7, mbC, sC, 0.5, 0, 1): ok = Not IsError(ph): chk = emptyChk
    If ok Then chk = EJC_CHECK(-98.0041250683277, ph, 77.7, mbC, sC, 0.5, 0, 1): ok = (Abs(ph - 89.95) < 0.01 And AcceptText(chk, 0.0000000001, 0.00000001) = "ADMISSIBLE")
    RecordTest ws, rw, "T4 valid HB point near 89.95 deg", "89.95, ADMISSIBLE", ErrText(ph), chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T5
    ph = EJC_PHI(400, 77.7, mbC, sC, 0.5, -1, 1 / 0.271): ok = IsErrCode(ph, xlErrNA)
    RecordTest ws, rw, "T5 closure 1, Chile, nu = 0.271", "#N/A (no root)", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T6
    ph = EJC_PHI(400, 0, mbC, sC, 0.5, 0, 1): ok = IsErrCode(ph, xlErrValue)
    RecordTest ws, rw, "T6 invalid input sigma_ci = 0", "#VALUE!", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T7
    hit = False: chk = emptyChk
    For i = 1 To 4
        ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, 0.0000000001, i)
        If Not IsError(ph) Then
            If Abs(ph - 39.89156042767222) < 0.000001 Then chk = EJC_CHECK(T7SN, ph, 120, 3.5, 0.002, 0.505, -1, 20): hit = (AcceptText(chk, 0.0000000001, 0.00000001) = "ADMISSIBLE")
        End If
    Next i
    RecordTest ws, rw, "T7 isolated eq. 4.21 counterexample (13.265 MPa)", "an ADMISSIBLE root at 39.8915604", "found: " & hit, chk, hit: rw = rw + 1: If hit Then nPass = nPass + 1
    ' T8, T9 (force-used state of T1)
    chk = EJC_CHECK(400, phT1, 77.7, mbC, sC, 0.5, 0, 1)
    If IsError(chk(5)) Then
        cV = 0
    Else
        cV = 1000 * ((chk(5) - chk(6)) / 2 * Cos(phT1 * PI_ / 180) - 0.4 * Tan(phT1 * PI_ / 180))
    End If
    chk = EJC_CHECK_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1): lbl = AcceptText(chk, 0.0000000001, 0.00000001): ok = (lbl = "ADMISSIBLE")
    RecordTest ws, rw, "T8 force-used state of T1", "ADMISSIBLE", lbl, chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    chk = EJC_CHECK_STATE(400, cV + 5, phT1, 77.7, mbC, sC, 0.5, 0, 1): lbl = AcceptText(chk, 0.0000000001, 0.00000001): ok = (lbl <> "ADMISSIBLE")
    RecordTest ws, rw, "T9 perturbed force-used state (c + 5 kPa)", "any REJECTED status", lbl, chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T10-T14 boundary handling
    ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, 0.0000000001, 0): ok = IsErrCode(ph, xlErrValue)
    RecordTest ws, rw, "T10 root index 0", "#VALUE!", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, 0.0000000001, -1): ok = IsErrCode(ph, xlErrValue)
    RecordTest ws, rw, "T11 root index -1", "#VALUE!", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, 0.0000000001, 2): ok = Not IsError(ph): chk = emptyChk
    If ok Then chk = EJC_CHECK(T7SN, ph, 120, 3.5, 0.002, 0.505, -1, 20): ok = (Abs(ph - 60.4228863398) < 0.000001 And AcceptText(chk, 0.0000000001, 0.00000001) <> "ADMISSIBLE")
    RecordTest ws, rw, "T12 root index 2 at T7 inputs", "60.4228863398, a REJECTED status", ErrText(ph), chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, 0, 1): ok = IsErrCode(ph, xlErrValue)
    RecordTest ws, rw, "T13 tolerance 0", "#VALUE!", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ph = EJC_PHI(T7SN, 120, 3.5, 0.002, 0.505, -1, 20, -0.0000000001, 1): ok = IsErrCode(ph, xlErrValue)
    RecordTest ws, rw, "T14 tolerance negative", "#VALUE!", ErrText(ph), emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T15-T17 branch / domain policy
    chk = EJC_CHECK_STATE(29819.0191175968, 25010.6868294027, -0.012480264998, 120, 3.5, 0.002, 0.505, -1, 20): lbl = AcceptText(chk, 0.0000000001, 0.00000001)
    ok = (lbl = "REJECTED: DOMAIN" Or lbl = "REJECTED: BRANCH")
    RecordTest ws, rw, "T15 negative-friction branch state (QA example)", "REJECTED: DOMAIN or BRANCH", lbl, chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    chk = EJC_CHECK(T7SN, -45, 120, 3.5, 0.002, 0.505, -1, 20): lbl = AcceptText(chk, 0.0000000001, 0.00000001): ok = (lbl = "REJECTED: DOMAIN")
    RecordTest ws, rw, "T16 candidate at -45 deg", "REJECTED: DOMAIN (controlled)", lbl, chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    chk = EJC_CHECK(T7SN, 90, 120, 3.5, 0.002, 0.505, -1, 20): lbl = AcceptText(chk, 0.0000000001, 0.00000001): ok = (lbl = "REJECTED: DOMAIN")
    RecordTest ws, rw, "T17 candidate at 90 deg", "REJECTED: DOMAIN (finite-angle policy)", lbl, chk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T18-T21: invalid acceptance tolerances must never certify a state
    chk = EJC_CHECK_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1)
    lbl = EJC_ACCEPT_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1, 0.0000000001, 0): ok = (lbl = "REJECTED: INPUT")
    RecordTest ws, rw, "T18 acceptance EJC tolerance zero", "REJECTED: INPUT", lbl, emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    lbl = EJC_ACCEPT_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1, 0.0000000001, -0.00000001): ok = (lbl = "REJECTED: INPUT")
    RecordTest ws, rw, "T19 acceptance EJC tolerance negative", "REJECTED: INPUT", lbl, emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    lbl = EJC_ACCEPT_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1, 0, 0.00000001): ok = (lbl = "REJECTED: INPUT")
    RecordTest ws, rw, "T20 acceptance angular tolerance zero", "REJECTED: INPUT", lbl, emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    lbl = EJC_ACCEPT_STATE(400, cV, phT1, 77.7, mbC, sC, 0.5, 0, 1, -0.0000000001, 0.00000001): ok = (lbl = "REJECTED: INPUT")
    RecordTest ws, rw, "T21 acceptance angular tolerance negative", "REJECTED: INPUT", lbl, emptyChk, ok: rw = rw + 1: If ok Then nPass = nPass + 1
    ' T22 is placed below the stable summary address so existing QA consumers keep reading G26.
    ph = EJC_PHI(-49.768342, 77.7, mbC, sC, 0.509923, 0, 1): ok = Not IsError(ph): chk = emptyChk
    If ok Then
        chk = EJC_CHECK(-49.768342, ph, 77.7, mbC, sC, 0.509923, 0, 1)
        ok = (Abs(ph - 71.27787820469955) < 0.0000001 And AcceptText(chk, 0.0000000001, 0.00000001) = "ADMISSIBLE")
    End If
    RecordTest ws, 28, "T22 wet a=.509923 tensile-apex scan regression", "71.2778782047, ADMISSIBLE", ErrText(ph), chk, ok
    If ok Then nPass = nPass + 1
    ws.Cells(26, 1).Value = "SUMMARY": ws.Cells(26, 3).Value = nPass & " of 22 passed"
    ws.Cells(26, 7).Value = IIf(nPass = 22, "ALL PASS", "FAILURES - do not use")
    ws.Columns("A:G").AutoFit
End Sub
