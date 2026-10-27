Attribute VB_Name = "Thesis_Final_QA"
Option Explicit
' GENG5512 thesis package v14.
' Native Excel finalization companion. Import alongside EJC_Eq421.bas v6.
' This tests the formula engine; it does not replace it with a primary VBA Hungr solver.
' v14a: QA_Errors now reports where error cells are; Thesis_ApplyPatch_v14a fixes the Excel-only
'       Geometry tangency error (one grid point at y = yc, x = xc - R). Run the patch once, then the QA.
Private Const PI_QA As Double = 3.141592653589793
Private Const TOL_F_QA As Double = 0.00000001
Private Const TOL_EQ_QA As Double = 0.000000001

Private Function QA_Number(ByVal valueIn As Variant) As Boolean
    QA_Number = False
    If IsError(valueIn) Then Exit Function
    If IsEmpty(valueIn) Then Exit Function
    If VarType(valueIn) = vbString Then Exit Function
    QA_Number = IsNumeric(valueIn)
End Function

Private Function QA_Near(ByVal valueIn As Variant, ByVal expected As Double) As Boolean
    If Not QA_Number(valueIn) Then Exit Function
    QA_Near = (Abs(CDbl(valueIn) - expected) <= TOL_F_QA * Application.Max(1#, Abs(expected)))
End Function

Private Function QA_Errors(Optional ByRef whereText As String) As Double
    ' Counts every formula or constant cell holding an Excel error value, in every worksheet,
    ' and lists the first eight blocks (sheet and address) in whereText.
    Dim sheetIn As Worksheet, bad As Range, kind As Variant, block As Range, nBlocks As Long
    whereText = ""
    For Each sheetIn In ThisWorkbook.Worksheets
        For Each kind In Array(xlCellTypeFormulas, xlCellTypeConstants)
            Set bad = Nothing
            On Error Resume Next
            Set bad = sheetIn.Cells.SpecialCells(kind, xlErrors)
            On Error GoTo 0
            If Not bad Is Nothing Then
                QA_Errors = QA_Errors + bad.CountLarge
                For Each block In bad.Areas
                    nBlocks = nBlocks + 1
                    If nBlocks <= 8 Then
                        If Len(whereText) > 0 Then whereText = whereText & "; "
                        whereText = whereText & sheetIn.Name & "!" & block.Address(False, False)
                    End If
                Next block
            End If
        Next kind
    Next sheetIn
    If nBlocks > 8 Then whereText = whereText & "; and " & (nBlocks - 8) & " more blocks"
End Function

Public Sub Thesis_ApplyPatch_v14a()
    ' One-time v14a workbook patch. Run once, check the message, then save the XLSM.
    ' Geometry row 61 is the last grid point on the y = yc line, at x = xc - R, which lies on the sphere
    ' by construction (radicand exactly 0). The x grid is built by repeated subtraction, so Excel can see a
    ' radicand of about -3E-10 there (published circle), SQRT fails and the point becomes "No real solution";
    ' the unguarded dip-angle cell E61 then shows #VALUE!. LibreOffice rounds that subtraction to 0, which
    ' is why the error appears only in Excel. The patch treats a rounding-level radicand
    ' (|radicand| <= 1E-9 R^2) as 0 in C61:D61 and wraps E42:E61 in IFERROR, as rows 63-481 already are.
    ' These cells feed the plots only; no factor-of-safety cell reads them.
    Dim g As Worksheet, radText As String, safeText As String, rowN As Long, okAll As Boolean
    Set g = ThisWorkbook.Worksheets("Geometry")
    radText = "$B$7^2-(B61-$B$5)^2-($A$41-$D$5)^2"
    safeText = "IF(ABS(" & radText & ")<=0.000000001*$B$7^2,0," & radText & ")"
    g.Range("E42:E61").Formula = "=IFERROR(DEGREES(ATAN((D42-D41)/(B14-B13))),""-"")"
    g.Range("C61").Formula = "=IFERROR(IF($C$5+SQRT(" & safeText & ")>($C$5+$B$7),""Outside sphere"",$C$5+SQRT(" & safeText & ")),""No real solution"")"
    g.Range("D61").Formula = "=IFERROR(IF($C$5-SQRT(" & safeText & ")>($C$5+$B$7),""Outside sphere"",$C$5-SQRT(" & safeText & ")),""No real solution"")"
    Application.CalculateFullRebuild
    okAll = True
    For rowN = 42 To 61
        If Left$(CStr(g.Range("E" & rowN).Formula), 9) <> "=IFERROR(" Then okAll = False
        If InStr(1, CStr(g.Range("E" & rowN).Formula), "B" & (rowN - 28) & "-B" & (rowN - 29), vbBinaryCompare) = 0 Then okAll = False
        If IsError(g.Range("E" & rowN).Value2) Then okAll = False
    Next rowN
    If IsError(g.Range("C61").Value2) Or IsError(g.Range("D61").Value2) Then okAll = False
    If okAll Then
        MsgBox "v14a patch applied to Geometry E42:E61, C61 and D61." & vbCrLf & _
            "Now run Thesis_CompletePackageQA, then save the XLSM.", vbInformation
    Else
        MsgBox "v14a patch did not verify. Do not save; close without saving and report this message.", vbExclamation
    End If
End Sub

Private Sub QA_Baseline(ByVal inputs As Worksheet)
    inputs.Range("C6").Value2 = 355.742311546406
    inputs.Range("E6").Value2 = 644.419294864537
    inputs.Range("C9").Value2 = 725.266997518736
    inputs.Range("D14").Value2 = 360
    inputs.Range("D16").Value2 = 1
    inputs.Range("E16").Value2 = 1 / Tan(50 * PI_QA / 180)
    inputs.Range("D20").Value2 = 7.2
    inputs.Range("D21").Value2 = 30
    inputs.Range("D22").Value2 = 27
    inputs.Range("D23:D25").Value2 = 0
    inputs.Range("P4").Value2 = 1
    inputs.Range("P5").Value2 = 1
    inputs.Range("P6").Value2 = 0.271
    inputs.Range("P7").Value2 = 77.7
    inputs.Range("P8").Value2 = 10
    inputs.Range("P9").Value2 = 42
    inputs.Range("P10").Value2 = 0
    inputs.Range("P14").Value2 = 0.5
    inputs.Range("P16").Value2 = 0.0000000001
    inputs.Range("P17").Value2 = 0.00000001
    inputs.Range("P18").Value2 = 0.00000001
    inputs.Range("P19").Value2 = 1
End Sub

Private Function QA_State(ByVal normal As Double, ByVal cohesion As Double, ByVal angleDeg As Double, _
    ByVal weight As Double, ByVal pressure As Double, ByVal area As Double, ByVal gVert As Double, _
    ByVal kVert As Double, ByVal fos As Double, ByVal useHB As Boolean, ByVal sciV As Double, _
    ByVal mbV As Double, ByVal sV As Double, ByVal aV As Double, ByVal pC As Double, ByVal qC As Double, _
    ByRef maxAngle As Double, ByRef maxFailure As Double, ByRef maxEq As Double, ByRef shear As Double) As Boolean
    Dim phi As Double, relative As Double, checks As Variant, label As String
    If weight <= 0 Or area <= 0 Or gVert <= 0 Or fos <= 0 Then Exit Function
    phi = angleDeg * PI_QA / 180
    If gVert + kVert * Tan(phi) / fos <= 0 Then Exit Function
    shear = cohesion + normal * Tan(phi)
    relative = (weight - ((normal + pressure) * area * gVert + shear * area * kVert / fos)) / weight
    If Abs(relative) > maxEq Then maxEq = Abs(relative)
    If Abs(relative) > TOL_EQ_QA Then Exit Function
    If useHB Then
        checks = EJC_CHECK_STATE(normal, cohesion, angleDeg, sciV, mbV, sV, aV, pC, qC)
        label = EJC_ACCEPT_STATE(normal, cohesion, angleDeg, sciV, mbV, sV, aV, pC, qC)
        If label <> "ADMISSIBLE" Then Exit Function
        If Abs(checks(1)) > maxAngle Then maxAngle = Abs(checks(1))
        If Abs(checks(2)) > maxFailure Then maxFailure = Abs(checks(2))
    End If
    QA_State = True
End Function

Private Function QA_ActualStates(ByRef maxAngle As Double, ByRef maxFailure As Double, ByRef maxEq As Double) As Boolean
    Dim inputs As Worksheet, columns As Worksheet, slices As Worksheet, rN As Long, cN As Long
    Dim sciV As Double, mbV As Double, sV As Double, aV As Double, pC As Double, qC As Double, useHB As Boolean
    Dim fos3 As Double, fos2 As Double, radius As Double, centreX As Double, kVert As Double, yV As Double
    Dim valuesIn As Variant, v As Variant, shear As Double, resistance As Double, driving As Double, countN As Long
    Set inputs = ThisWorkbook.Worksheets("Interface")
    Set columns = ThisWorkbook.Worksheets("Iteration 34")
    Set slices = ThisWorkbook.Worksheets("2D Bishop (y = 0)")
    useHB = (inputs.Range("P4").Value2 = 1)
    sciV = inputs.Range("P7").Value2: mbV = inputs.Range("P11").Value2
    sV = inputs.Range("P12").Value2: aV = inputs.Range("P15").Value2
    Select Case inputs.Range("P5").Value2
        Case 1: pC = -1: qC = 1 / inputs.Range("P6").Value2
        Case 2: pC = inputs.Range("P6").Value2: qC = pC
        Case 3: pC = 0: qC = 1
        Case Else: Exit Function
    End Select
    If Not QA_Number(inputs.Range("M34").Value2) Then Exit Function
    If Not QA_Number(slices.Range("PF13").Value2) Then Exit Function
    fos3 = inputs.Range("M34").Value2: fos2 = slices.Range("PF13").Value2
    radius = inputs.Range("C9").Value2: centreX = inputs.Range("C6").Value2
    If fos3 <= 0 Or fos2 <= 0 Or radius <= 0 Then Exit Function
    resistance = 0: driving = 0: countN = 0
    For cN = 4 To 24
        For rN = 11 To 31
            If QA_Number(columns.Cells(rN + 77, cN).Value2) Then
                If columns.Cells(rN + 77, cN).Value2 > 0 Then
                    valuesIn = Array(columns.Cells(rN + 706, cN).Value2, columns.Cells(rN + 238, cN).Value2, _
                        columns.Cells(rN + 212, cN).Value2, columns.Cells(rN + 77, cN).Value2, _
                        columns.Cells(rN + 846, cN).Value2, columns.Cells(rN, cN).Value2, columns.Cells(rN + 25, cN).Value2)
                    For Each v In valuesIn
                        If Not QA_Number(v) Then Exit Function
                    Next v
                    yV = columns.Cells(9, cN).Value2
                    If radius * radius - yV * yV <= 0 Then Exit Function
                    kVert = (centreX - columns.Cells(rN, 3).Value2) / Sqr(radius * radius - yV * yV)
                    If Not QA_State(CDbl(valuesIn(0)), CDbl(valuesIn(1)), CDbl(valuesIn(2)), CDbl(valuesIn(3)), _
                        CDbl(valuesIn(4)), CDbl(valuesIn(5)), CDbl(valuesIn(6)), kVert, fos3, useHB, sciV, mbV, sV, aV, pC, qC, _
                        maxAngle, maxFailure, maxEq, shear) Then Exit Function
                    countN = countN + 1
                    resistance = resistance + shear * CDbl(valuesIn(5))
                    driving = driving + CDbl(valuesIn(3)) * kVert
                End If
            End If
        Next rN
    Next cN
    If countN <> inputs.Range("P22").Value2 Or countN = 0 Or driving <= 0 Then Exit Function
    If Not QA_Near(resistance / driving, fos3) Then Exit Function
    resistance = 0: driving = 0: countN = 0
    For rN = 14 To 53
        If slices.Cells(rN, 6).Value2 = 1 Then
            valuesIn = Array(slices.Range("PH" & rN).Value2, slices.Range("PD" & rN).Value2, slices.Range("PC" & rN).Value2, _
                slices.Range("G" & rN).Value2, slices.Range("L" & rN).Value2, slices.Range("J" & rN).Value2, _
                slices.Range("H" & rN).Value2, slices.Range("I" & rN).Value2)
            For Each v In valuesIn
                If Not QA_Number(v) Then Exit Function
            Next v
            If Not QA_State(CDbl(valuesIn(0)), CDbl(valuesIn(1)), CDbl(valuesIn(2)), CDbl(valuesIn(3)), _
                CDbl(valuesIn(4)), CDbl(valuesIn(5)), CDbl(valuesIn(6)), CDbl(valuesIn(7)), fos2, useHB, sciV, mbV, sV, aV, pC, qC, _
                maxAngle, maxFailure, maxEq, shear) Then Exit Function
            countN = countN + 1
            resistance = resistance + shear * CDbl(valuesIn(5))
            driving = driving + CDbl(valuesIn(3)) * CDbl(valuesIn(7))
        End If
    Next rN
    If countN <> slices.Range("B4").Value2 Or countN = 0 Or driving <= 0 Then Exit Function
    If Not QA_Near(resistance / driving, fos2) Then Exit Function
    QA_ActualStates = True
End Function

Public Sub Thesis_RunFinalQA()
    Dim inputs As Worksheet, result As Worksheet, slices As Worksheet, columns As Worksheet
    Dim addresses As Variant, stored() As Variant, idxN As Long, testN As Long, rowN As Long, nPass As Long
    Dim oldCalc As XlCalculation, oldEvents As Boolean, oldScreen As Boolean, originalFinal As Variant, perturbActive As Boolean
    Dim expected3 As Double, expected2 As Double, nameTest As String, ok As Boolean, stateOK As Boolean, unitsOK As Boolean
    Dim maxAngle As Double, maxFailure As Double, maxEq As Double, errorsN As Double, approval As String, repair As String
    Dim savedInputs As Boolean, runtimeText As String, restoredOK As Boolean
    Dim crossOK As Boolean, phCross As Variant, errWhere As String
    On Error GoTo Failed
    Set inputs = ThisWorkbook.Worksheets("Interface")
    Set slices = ThisWorkbook.Worksheets("2D Bishop (y = 0)")
    Set columns = ThisWorkbook.Worksheets("Iteration 34")
    addresses = Array("C6", "E6", "C9", "D14", "D16", "E16", "D20", "D21", "D22", "D23", "D24", "D25", _
        "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P14", "P16", "P17", "P18", "P19")
    ReDim stored(LBound(addresses) To UBound(addresses))
    For idxN = LBound(addresses) To UBound(addresses): stored(idxN) = inputs.Range(addresses(idxN)).Formula: Next idxN
    savedInputs = True
    oldCalc = Application.Calculation: oldEvents = Application.EnableEvents: oldScreen = Application.ScreenUpdating
    Application.Calculation = xlCalculationManual: Application.EnableEvents = False: Application.ScreenUpdating = False
    On Error Resume Next
    Set result = ThisWorkbook.Worksheets("Native Workbook QA")
    On Error GoTo Failed
    If result Is Nothing Then
        Set result = ThisWorkbook.Worksheets.Add(After:=ThisWorkbook.Worksheets(ThisWorkbook.Worksheets.Count))
        result.Name = "Native Workbook QA"
    Else
        approval = CStr(result.Range("B18").Value2): repair = CStr(result.Range("B23").Value2)
    End If
    result.Cells.Clear
    result.Range("A1").Value2 = "v14 native Excel workbook QA - Excel " & Application.Version & " - " & Format(Now, "yyyy-mm-dd hh:nn:ss")
    result.Range("A2").Value2 = "Workbook: " & ThisWorkbook.Name & "; v14a / EJC_Eq421 revision 6; F 1E-8 relative, angle 1E-10 rad, EJC 1E-8 MPa, equilibrium 1E-9 relative."
    result.Range("A3:M3").Value = Array("Test", "3D expected", "3D actual/internal", "2D expected", "2D actual/internal", _
        "3D status", "2D status", "Error cells", "Max angle rad", "Max EJC MPa", "Max local eq/W", "Result", "Notes")
    EJC_RunValidation
    unitsOK = (ThisWorkbook.Worksheets("VBA Validation").Range("G26").Value2 = "ALL PASS" And _
        ThisWorkbook.Worksheets("VBA Validation").Range("C26").Value2 = "22 of 22 passed")
    nPass = 0
    For testN = 1 To 9
        Application.StatusBar = "Thesis final QA: case " & testN & " of 9 (full recalculation)."
        DoEvents
        QA_Baseline inputs
        expected3 = 0: expected2 = 0: maxAngle = 0: maxFailure = 0: maxEq = 0: stateOK = False
        Select Case testN
            Case 1: nameTest = "Default thesis closure: rejection"
            Case 2: nameTest = "HB dry a=0.5": inputs.Range("P5").Value2 = 3
                expected3 = 2.5536684200376: expected2 = 2.31025502757599
            Case 3: nameTest = "HB wet a=0.5": inputs.Range("P5").Value2 = 3
                inputs.Range("D24").Value2 = 50: inputs.Range("D25").Value2 = 360
                expected3 = 1.20012184467095: expected2 = 1.03587083245705
            Case 4: nameTest = "HB dry a=0.509923": inputs.Range("P5").Value2 = 3: inputs.Range("P14").Value2 = 0.509923
                expected3 = 2.47603808311863: expected2 = 2.24206326532995
            Case 5, 9: nameTest = "HB wet a=0.509923": inputs.Range("P5").Value2 = 3: inputs.Range("P14").Value2 = 0.509923
                inputs.Range("D24").Value2 = 50: inputs.Range("D25").Value2 = 360
                expected3 = 1.15842745018062: expected2 = 1.00046615780691
            Case 6: nameTest = "Published circle, conventional HB": inputs.Range("P5").Value2 = 3
                inputs.Range("C6").Value2 = 207.28: inputs.Range("E6").Value2 = 586.53: inputs.Range("C9").Value2 = 622.08
                expected3 = 2.14193739499583: expected2 = 2.00387341089056
            Case 7: nameTest = "MC zero-friction numerical check": inputs.Range("P4").Value2 = 0
                inputs.Range("D20").Value2 = 100: inputs.Range("D21").Value2 = 0
                expected3 = 0.255354533635107: expected2 = 0.190263799071471
            Case 8: nameTest = "Moscarda wet diagnostic: result withheld": inputs.Range("P5").Value2 = 2: inputs.Range("P19").Value2 = 0
                inputs.Range("P14").Value2 = 0.509923: inputs.Range("D24").Value2 = 50: inputs.Range("D25").Value2 = 360
                expected3 = 1.04327747791523: expected2 = 0.90176054371878
        End Select
        Application.CalculateFullRebuild
        If testN = 2 Then
            phCross = EJC_PHI(columns.Range("E182").Value2, inputs.Range("P7").Value2, inputs.Range("P11").Value2, _
                inputs.Range("P12").Value2, inputs.Range("P15").Value2, 0, 1)
            crossOK = False
            If Not IsError(phCross) Then crossOK = (Abs(CDbl(phCross) - columns.Range("E234").Value2) <= 0.00000001)
            crossOK = crossOK And (EJC_ACCEPT_STATE(columns.Range("E728").Value2, columns.Range("E260").Value2, _
                columns.Range("E234").Value2, inputs.Range("P7").Value2, inputs.Range("P11").Value2, _
                inputs.Range("P12").Value2, inputs.Range("P15").Value2, 0, 1) = "ADMISSIBLE")
        End If
        If testN = 9 Then
            nameTest = "Final stress +1 kPa: result withheld"
            originalFinal = columns.Range("J730").Formula: perturbActive = True
            columns.Range("J730").Value2 = CDbl(columns.Range("J730").Value2) + 1
            Application.CalculateFullRebuild
        End If
        errorsN = QA_Errors(errWhere)
        If testN = 1 Then
            ok = (inputs.Range("P31").Value2 = "NOT ACCEPTED" And inputs.Range("D31").Value2 = "n.a." And _
                inputs.Range("P33").Value2 = "n.a." And inputs.Range("P22").Value2 = 76 And inputs.Range("P23").Value2 = 0)
        ElseIf testN = 8 Then
            ok = (Left(CStr(inputs.Range("P31").Value2), 15) = "DIAGNOSTIC ONLY" And inputs.Range("D31").Value2 = "n.a." And _
                inputs.Range("P33").Value2 = "n.a." And inputs.Range("P28").Value2 > inputs.Range("P17").Value2)
            ok = ok And QA_Near(inputs.Range("M34").Value2, expected3) And QA_Near(slices.Range("PF13").Value2, expected2)
        ElseIf testN = 9 Then
            ok = (inputs.Range("P31").Value2 = "NOT ACCEPTED" And inputs.Range("D31").Value2 = "n.a." And inputs.Range("P26").Value2 > 0)
        Else
            stateOK = QA_ActualStates(maxAngle, maxFailure, maxEq)
            ok = (inputs.Range("P31").Value2 = "ACCEPTED" And slices.Range("E6").Value2 = "ACCEPTED")
            ok = ok And QA_Near(inputs.Range("D31").Value2, expected3) And QA_Near(inputs.Range("P33").Value2, expected2) And stateOK
        End If
        ok = ok And (errorsN = 0)
        rowN = testN + 3
        result.Cells(rowN, 1).Value2 = nameTest
        If expected3 > 0 Then result.Cells(rowN, 2).Value2 = expected3 Else result.Cells(rowN, 2).Value2 = "withheld"
        result.Cells(rowN, 3).Value2 = inputs.Range("M34").Value2
        If expected2 > 0 Then result.Cells(rowN, 4).Value2 = expected2 Else result.Cells(rowN, 4).Value2 = "withheld"
        result.Cells(rowN, 5).Value2 = slices.Range("PF13").Value2
        result.Cells(rowN, 6).Value2 = inputs.Range("P31").Value2: result.Cells(rowN, 7).Value2 = slices.Range("E6").Value2
        result.Cells(rowN, 8).Value2 = errorsN
        If stateOK Then
            If testN = 7 Then
                result.Cells(rowN, 9).Value2 = "n.a. (MC)": result.Cells(rowN, 10).Value2 = "n.a. (MC)"
            Else
                result.Cells(rowN, 9).Value2 = maxAngle: result.Cells(rowN, 10).Value2 = maxFailure
            End If
            result.Cells(rowN, 11).Value2 = maxEq
        Else
            result.Range("I" & rowN & ":K" & rowN).Value2 = "not evaluated as an accepted state"
        End If
        result.Cells(rowN, 12).Value2 = IIf(ok, "PASS", "FAIL")
        If testN = 7 Then result.Cells(rowN, 13).Value2 = "Robustness test; no calibrated field material."
        If testN = 9 Then
            result.Cells(rowN, 13).Value2 = "One actual final stress perturbed, then restored."
            columns.Range("J730").Formula = originalFinal: perturbActive = False
        End If
        If errorsN > 0 Then result.Cells(rowN, 13).Value2 = Trim$(CStr(result.Cells(rowN, 13).Value2) & " Error cells at: " & errWhere)
        If ok Then nPass = nPass + 1
    Next testN
    result.Range("A13").Value2 = "VBA solver vs worksheet phi (dry-HB E22)"
    result.Range("L13").Value2 = IIf(crossOK, "PASS", "FAIL")
    If crossOK Then nPass = nPass + 1
    result.Range("A15").Value2 = "NUMERICAL SUMMARY"
    result.Range("C15").Value2 = nPass & " of 10 workbook checks; unit fixtures: " & IIf(unitsOK, "22 of 22 PASS", "FAIL")
    result.Range("L15").Value2 = IIf(nPass = 10 And unitsOK, "ALL PASS", "FAILURES - do not submit as verified")
    result.Range("A18").Value2 = "Supervisor accepts documented C1/software role and negative Chile outcome?"
    result.Range("B18").Value2 = IIf(approval = "APPROVED", approval, "OPEN")
    result.Range("A19").Value2 = "Implemented software role": result.Range("B19").Value2 = "Primary VBA assembly supplied in Thesis_Columns; worksheets retained as independent comparison."
    result.Range("A20").Value2 = "C1 / stress closure": result.Range("B20").Value2 = "Eq. 4.21 default retained; sigma2 = sigma3/nu - sigma1."
    result.Range("A21").Value2 = "Chile thesis-closure outcome": result.Range("B21").Value2 = "No admissible FoS; no fallback strength."
    result.Range("A22").Value2 = "Worksheet QA gate (overall gate: Completion Status)"
    result.Range("B22").Formula = "=IF(AND(L15=""ALL PASS"",B18=""APPROVED"",B23=""YES""),""WORKSHEET QA GATES RECORDED - inspect Completion Status"",""OPEN - native QA/method agreement/repair check required"")"
    result.Range("A23").Value2 = "Manually confirmed: opened with no Excel repair prompt?"
    result.Range("B23").Value2 = IIf(repair = "YES", repair, "OPEN")
    result.Range("A25").Value2 = "Primary VBA assembly: Thesis_Columns. This sheet tests the worksheet comparison. Native Engine QA tests the primary assembly; Completion Status combines both records."
    result.Columns("A:M").AutoFit
    GoTo Restore
Failed:
    runtimeText = "Error " & Err.Number & ": " & Err.Description
    On Error Resume Next
    If Not result Is Nothing Then
        result.Range("L15").Value2 = "FAIL - interrupted"
        result.Range("A16").Value2 = runtimeText
    End If
Restore:
    On Error Resume Next
    restoredOK = True
    If perturbActive Then
        Err.Clear: columns.Range("J730").Formula = originalFinal
        If Err.Number <> 0 Then restoredOK = False
    End If
    If savedInputs Then
        For idxN = LBound(addresses) To UBound(addresses)
            Err.Clear: inputs.Range(addresses(idxN)).Formula = stored(idxN)
            If Err.Number <> 0 Then restoredOK = False
            Err.Clear
            If CStr(inputs.Range(addresses(idxN)).Formula) <> CStr(stored(idxN)) Then restoredOK = False
            If Err.Number <> 0 Then restoredOK = False
        Next idxN
        Err.Clear
        Application.CalculateFullRebuild
        If Err.Number <> 0 Then restoredOK = False
        Application.Calculation = oldCalc: Application.EnableEvents = oldEvents: Application.ScreenUpdating = oldScreen
    End If
    If Not result Is Nothing Then result.Range("A16").Value2 = "Input/formula restoration: " & IIf(restoredOK, "exact match", "FAILED")
    If Not restoredOK Then
        runtimeText = runtimeText & " Input restoration/recalculation failed; review the inputs before use."
        If Not result Is Nothing Then result.Range("L15").Value2 = "FAIL - restoration incomplete"
    End If
    Application.StatusBar = False
    On Error GoTo 0
    If Len(runtimeText) > 0 Then MsgBox runtimeText & vbCrLf & "Tests stopped; do not treat this record as verified.", vbExclamation
    If Len(runtimeText) = 0 Then MsgBox "Native records written. Save the XLSM after reviewing both validation sheets and the submission gate.", vbInformation
End Sub
