# -*- coding:UTF-8 -*-

# Этот скрипт предназначен для достраивания полугоризонталей. В качестве входных данных используются
# горизонтали, содержащие значения высоты

# ПЕРЕД ЗАПУСКОМ ОБРАТИТЕ ВНИМАНИЕ НА СТРОКУ 55

# подключение необходимых библиотек 
import arcpy
from arcpy.sa import *
import arcpy.cartography as CA



def contour_prolongation(input_contours, 
                         height_field, 
                         intermediate_raster_cellsize, 
                         output_contours):
    # Описание свойств слоя с горизонталями
    desc = arcpy.Describe(input_contours)
    xmin = desc.extent.XMin
    xmax = desc.extent.XMax
    ymin = desc.extent.YMin
    ymax = desc.extent.YMax
    # экстент слоя
    arcpy.env.extent = desc.extent

    # экстент по Х, по У и полный
    Xext = xmax - xmin
    Yext = ymax - ymin
    full_ext = round((Xext + Yext), 0)

    # размер ячейки 
    cellSize = 3

    arcpy.AddMessage("Creation of the rasters with constant value")
    # Создаем 2 растра постоянного значения (constant raster) для записи значений
    # Задаем очень большие (вряд ли достижимые) значения евклидового расстояния. Экстент по Х+ экстент по У
    minVal_1 = CreateConstantRaster(full_ext, "FLOAT", cellSize) 
    minVal_1.save("minVal_1")

    # перевод в массив
    # здесь будет записано первое минимальное значение после цикла
    minVal_1_numpy_array = arcpy.RasterToNumPyArray(minVal_1) 

    # то же для второго растра
    minVal_2 = CreateConstantRaster(full_ext, "FLOAT", cellSize) 
    minVal_2.save("minVal_2")

    # здесь будет записано второе минимальное значение после цикла
    minVal_2_numpy_array = arcpy.RasterToNumPyArray(minVal_2) 

    arcpy.AddMessage("Calculation of the contour interval")
    # НАЗВАНИЕ ПОЛЯ, СОДЕРЖАЩЕГО ВЫСОТЫ ГОРИЗОНТАЛЕЙ
    field = "Habs"
    # ЕСЛИ ПОЛЕ С ВЫСОТАМИ НАЗЫВАЕТСЯ НЕ 'Habs', ТО РАСКОММЕНТИРУЙТЕ СЛЕДУЮЩУЮ СТРОКУ
    # arcpy.AlterField_management(fc, field, 'Habs') # переименование поля с высотами

    # поиск высот горизонталей
    Hvalues = [row[0] for row in arcpy.da.SearchCursor(input_contours, "Habs")]
    uniqueValues = set(Hvalues)
    # конвертация в список, уникальные значения высот
    Hvalues = list(uniqueValues) 

    heights = Hvalues
    heights_unique = list(set(heights))
    # сортировка уникальных значений высот по возрастанию
    heights_unique.sort()
    arcpy.AddMessage("Unique z-values:")
    print(heights_unique)

    # список, в который будут записана высота сечения между соседними горизонталями
    steps = []
    # проход по значениям высот и нахождение разницы между следующим и текущим 
    # значениями - это высота сечения
    for i in range(len(heights_unique) - 1):
        step = heights_unique[i + 1] - heights_unique[i]
        # добавление в список
        steps.append(step)
    # print(steps)
    steps_unique = list(set(steps))
    # сортировка значений по возрастанию
    steps_unique.sort()
    # последнее значение списка соответствует высоте сечения рельефа
    contour_step = steps_unique[-1]
    arcpy.AddMessage("The contour interval is {0} m".format(contour_step))

    arcpy.AddMessage("Creation of the new field to divide contour lines into intermediate and supplementary")
    arcpy.AddField_management(input_contours, "line_type", "TEXT", field_alias = "line_type")
    fields = ["Habs", "line_type"]
    with arcpy.da.UpdateCursor(input_contours, fields) as cursor:
        for row in cursor:
            # если остаток от деления высоты горизонтали на высоту сечения рельефа равен 0, то это
            # основная горизонталь
            if(row[0] % contour_step == 0):
                row[1] = "intermediate"
            else:
                # если остаток от деления высоты горизонтали на высоту сечения рельефа не равен 0, 
                # то это полугоризонталь
                row[1] = "supplementary"
            cursor.updateRow(row)

    # Создаем слой из исходного файла
    arcpy.MakeFeatureLayer_management(input_contours, "lyr") 

    # Делаем выборку только основных горизонталей, так найдем срединные линии (будущие полугоризонтали)
    arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", "\"line_type\" = 'intermediate'")


    arcpy.AddMessage("Euclidean distance calculation")
    # максимальное расстояние, больше которого ячейкам будет присвоено значение "no data"
    maxDistance = full_ext
    # евклидово расстояние нужно для создания срединных линий между 2 основными горизонталями
    outEucDistance = EucDistance("lyr", maxDistance, cellSize) 
    # сохранение 
    outEucDistance.save("EuDist")

    # слой подвыборки
    subset_layer = "lyr"
    # размеры растров
    nrows = outEucDistance.height
    ncols = outEucDistance.width

    # Список уникальных значений высот для подвыборки
    subHvalues = [row[0] for row in arcpy.da.SearchCursor(subset_layer, "Habs")]
    uniqueSubHValues = set(subHvalues)
    subHvalues = list(uniqueSubHValues) # конвертация в список, уникальные значения высот
    arcpy.AddMessage("Unique z-values:")
    print(subHvalues)

    arcpy.AddMessage("Comparison new euclidean distance with primary one")
    # В этом цикле будет происходить сравнение значений нового евклидового расстояния (ЕР) с изначальным
    for i in subHvalues:
        # Делаем выборку внутри выборки
        arcpy.SelectLayerByAttribute_management("lyr", "NEW_SELECTION", "Habs = %s" % (i)) # перебор по высотам; 
        # %s - формат строки, далее подстановка
        # новое ЕР
        outEucDistance = EucDistance("lyr", maxDistance, cellSize)  
        EucDist_numpy = arcpy.RasterToNumPyArray(outEucDistance)
        # сравнение нового ЕР с изначальным
        for row in range(nrows):
            for col in range(ncols):
                # если ЕР для данной точки меньше того, что в растре постоянного значения (1), 
                # то перезаписываем это значение в него. Затем переносим это значение в растр постоянного 
                # значения (2)
                if(EucDist_numpy[row, col] <= minVal_1_numpy_array[row, col]): # обращение к ячейке
                    minVal_2_numpy_array[row, col] = minVal_1_numpy_array[row, col]
                    minVal_1_numpy_array[row, col] = EucDist_numpy[row, col] 
                # Если ЕР для данной точки больше ЕР, записанного в растре постоянного значения (1), 
                # и меньше ЕР, записанного в растре постоянного значения (2), то новое значение 
                # пишем в растр постоянного значения (2)
                elif(EucDist_numpy[row, col] > minVal_1_numpy_array[row, col] and 
                EucDist_numpy[row, col] < minVal_2_numpy_array[row, col]):
                    minVal_2_numpy_array[row, col] = EucDist_numpy[row, col]

    arcpy.AddMessage("Middle lines creation")                
    # перезаписанные растры постоянного значения с учетом евклидового расстояния
    # перевод из массива в растровый формат
    minVal_1_new = arcpy.NumPyArrayToRaster(minVal_1_numpy_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    minVal_1_new.save("minVal_1_new")

    minVal_2_new = arcpy.NumPyArrayToRaster(minVal_2_numpy_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    minVal_2_new.save("minVal_2_new")
            
    # Разница между 1 и 2 растрами минимальных значений: ячейки со значениями, близкими к 0, показывают 
    # среднюю линию между 2 основными горизонталями
    diff_array = minVal_1_numpy_array - minVal_2_numpy_array
    diff = arcpy.NumPyArrayToRaster(diff_array, arcpy.Point(xmin, ymin), cellSize, cellSize)
    # Преобразование в целочисленное значение
    diff_int = Int(diff)  

    # установка нулевых значений, чтобы отделить срединную линию от фона
    # значение подобрано эмпирически
    set0 = SetNull(diff_int, diff_int, '"VALUE" < -5')
    # Преобразование в целочисленное значение
    set0_int = Int(set0)

    arcpy.AddMessage("Middle line processing") 
    # Утончение растра
    # фоновые значения - nodata, т.к. на предыдущем этапе они заданы; скругление углов (round)
    ThinRaster = Thin(set0_int, background_value = "NODATA", filter = "NO_FILTER", corners = "ROUND") 
    # Перевод растра в линии
    midlines = "midlines"
    # фоновые значения = 0, упрощать геометрию (SIMPLIFY)
    arcpy.RasterToPolyline_conversion(ThinRaster, midlines, background_value = "ZERO", 
    simplify = "SIMPLIFY", raster_field = "VALUE") 

    # Оптимизация формы линий
    # параметры подобраны эмпирически
    midlines_smpl = "midlines_smpl"
    CA.SimplifyLine(midlines, midlines_smpl, "BEND_SIMPLIFY", tolerance = 25) # диаметр сглаживания - 25 м
    midlines_smooth = "midlines_smooth"
    CA.SmoothLine(midlines_smpl, midlines_smooth, "PAEK", tolerance = 5)

    arcpy.AddMessage("Selection of middle lines") 
    arcpy.MakeFeatureLayer_management(input_contours, "suppl_lines")  
    # Отделение полугоризонталей 
    suppl_lines = arcpy.SelectLayerByAttribute_management("suppl_lines", "NEW_SELECTION",
    "\"line_type\" = 'supplementary'")

    arcpy.MakeFeatureLayer_management(midlines_smooth, "midlines") 
    # Выбор из срединных линий тех, которые пересекаются с исходными полугоризонталями
    midl_intersect = arcpy.SelectLayerByLocation_management("midlines", "INTERSECT", suppl_lines)

    arcpy.AddMessage("Intersect middle lines with supplementary lines") 
    # Пересечение исходных полугоризонталей со срединными линиями, на выходе получаем точки
    for_inters = [suppl_lines, midl_intersect]
    inters_point = "intersection_point"
    arcpy.Intersect_analysis(for_inters, inters_point, "ALL", output_type = "point")

    # Перевод точек пересечения в отдельные объекты
    inters_point_single = "intersection_point_single"
    arcpy.MultipartToSinglepart_management(inters_point, inters_point_single)

    # Разрезать срединные линии по точкам пересечений
    split_lns = "splitted_lines"
    arcpy.SplitLineAtPoint_management(midl_intersect, inters_point_single, split_lns, 
                                      search_radius= "1 Meters")

    arcpy.AddMessage("Spatial join of intersection points and middle lines")
    # Пространственное соединение точек пересечения и отобранных срединных линий, чтобы получить
    # поле высот (Habs)
    sp_join_midl_intersect = "sp_join_midl_intersect"
    arcpy.SpatialJoin_analysis(split_lns, inters_point_single, match_option= "INTERSECT", 
                               out_feature_class = sp_join_midl_intersect, join_operation = "JOIN_ONE_TO_ONE", 
                               join_type = "KEEP_ALL", search_radius= "1 Meters")

    arcpy.AddMessage("Processing of every segment of middle lines")
    # выборка необходимых слоев и полей
    arcpy.MakeFeatureLayer_management(sp_join_midl_intersect, "split_lns")
    arcpy.MakeFeatureLayer_management(inters_point_single, "intersect_pts")
    arcpy.da.SearchCursor(inters_point_single, "FID_new_lines_full")
    midlines_split = arcpy.da.SearchCursor(sp_join_midl_intersect, "OBJECTID")

    arcpy.AddMessage("Intersection point selection")
    for i in midlines_split:
        # Выделение каждого отрезка срединной линии
        selection = "OBJECTID = {}".format(i[0])       
        lns = arcpy.SelectLayerByAttribute_management("split_lns", "NEW_SELECTION", selection)
        # выбор точек, которые пересекаются с отрезками
        sel = arcpy.SelectLayerByLocation_management("intersect_pts", "INTERSECT", lns, 
        search_distance = "1 Meters", selection_type = "NEW_SELECTION")
        # подсчет количества выбранных точек
        counts = int(arcpy.GetCount_management(sel)[0])
        
        # Если точка 1, то сохраняем линию в новый слой. Если 2 точки, то получаем ID исходной 
        # полугоризонтали. Если ID у 2 точек разные, то линию сохраняем
        if(counts == 1):
            # проверка существования слоя. Если его нет, создаем. Если есть, то дополняем
            if arcpy.Exists("selected_lns"):
                arcpy.Append_management(lns, "selected_lns", "NO_TEST")
            else:
                arcpy.FeatureClassToFeatureClass_conversion(lns, arcpy.env.workspace, "selected_lns")
        if(counts == 2):
            # выбор тех срединных линий, у точек пересечения которых разные идентификаторы, т.е. разные
            # полугоризонтали, которые соединяем
            ids = [row[0] for row in arcpy.da.SearchCursor(sel, "OBJECTID")]
            if(ids[0] != ids[1]):
                if arcpy.Exists("selected_lns"):
                    arcpy.Append_management(lns, "selected_lns", "NO_TEST")
                else:
                    arcpy.FeatureClassToFeatureClass_conversion(lns, arcpy.env.workspace, "selected_lns")
                    


    arcpy.AddMessage("Creation of Voronoy diagram")
    # создание диаграмм Вороного для всех исходных горизонталей, чтобы далее сделать обрезку срединных линий
    # сначала необходимо уплотнить точки в изолиниях и перевести изолинии в точки
    # копирование слоя, т.к. при уплотнении новый слой не создается
    full_lns_copy = "full_lns_copy"
    arcpy.CopyFeatures_management(input_contours, full_lns_copy)
    arcpy.Densify_edit(full_lns_copy, "DISTANCE")

    pts_from_all_lines = "pts_from_all_lines"
    arcpy.FeatureVerticesToPoints_management(full_lns_copy, pts_from_all_lines, "ALL")
    
    # создание диаграмм Вороного
    thiess_poly = "thiess_poly"
    arcpy.CreateThiessenPolygons_analysis(pts_from_all_lines, thiess_poly, fields_to_copy = "ALL")
    arcpy.MakeFeatureLayer_management(thiess_poly, "thiess_poly")
    # Выбор полигонов, относящихся к полугоризонталям
    suppl_poly = arcpy.SelectLayerByAttribute_management("thiess_poly", "NEW_SELECTION",
    "\"line_type\" = 'supplementary'")
    # Обрезка срединных линий по выбранным полигонам
    clipped_lns_by_polyg = "clipped_lns_by_polyg"
    arcpy.Erase_analysis("selected_lns", suppl_poly, clipped_lns_by_polyg)
    arcpy.Delete_management(full_lns_copy)
    arcpy.Delete_management(pts_from_all_lines)
    arcpy.Delete_management(thiess_poly)

    arcpy.AddMessage("Merge middli lines with supplemenatry lines")
    # Соединение срединных линий с полугоризонталями
    full_lines_copy = 'full_lines_copy'
    arcpy.CopyFeatures_management(input_contours, full_lines_copy)
    fieldMappings = arcpy.FieldMappings()
    # добавление всех атрибутов двух слоев
    fieldMappings.addTable(full_lines_copy)
    fieldMappings.addTable(clipped_lns_by_polyg)
    fldMap = arcpy.FieldMap()
    fldMap.addInputField(full_lines_copy, "Habs")
    fieldMappings.addFieldMap(fldMap)

    # удаление всех атрибутов, кроме высоты и ID
    for field in fieldMappings.fields:
        if field.name not in ["OBJECTID", "Habs"]:
            fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex(field.name))

    joined_contours = "joined_contours"
    arcpy.Merge_management([full_lines_copy, clipped_lns_by_polyg], joined_contours, fieldMappings)
    # удаление дублирующегося поля с высотой
    arcpy.DeleteField_management("joined_contours", "Habs_1")

    # создание конечных точек у изолиний, чтобы их далее соединить
    end_pts_from_lns = "end_pts_from_lns"
    arcpy.FeatureVerticesToPoints_management(joined_contours, end_pts_from_lns, "BOTH_ENDS")
    # поиск ближайших концов изолиний внутри слоя
    # input и near - конечные точки
    arcpy.Near_analysis(end_pts_from_lns, end_pts_from_lns, search_radius = "500 Meters")

    # если в полученном слое OBJECTID_1[0] = NEAR_FID[1] и OBJECTID_1[1] = NEAR_FID[0], 
    # то соединяем эти точки
    arcpy.AddMessage("Find point pairs:")
    fields = ["OBJECTID", "NEAR_FID", "Habs", "NEAR_DIST"]
    fids_list = []
    # список с точками, которые не подходят
    skip_list = []
    for k in arcpy.da.SearchCursor("end_pts_from_lns", fields):
        # текущая строка, проверка неравенства расстояния нулю (строгое равенство может не сработать)
        if(k[0] in skip_list or k[3] < 1e-7):
            continue
        # выбор подходящей точки
        # условие: id ближайшей точки (NEAR_FID) = id текущей точки (OBJECTID_1). Если в обратную сторону
        # условие верно, то такие 2 точки считаются ближайшими и подлежат соединению
        where_clause = "OBJECTID = %s" % (k[1]) # подстановка
        for j in arcpy.da.SearchCursor("end_pts_from_lns", fields, where_clause):
            if(j[1] == k[0] and k[2] == j[2]):  # проверить равенство высот
                fids_list.append([k[0], j[0], k[2]])
                # для избежания дублирования пар точек вторую пару пропускаем
                skip_list.append(j[0])
    # print(fids_list)

    arcpy.AddMessage("Convert points to lines and write z-values")
    # копирование слоя, чтобы не испортить при присоединении фрагментов линий
    pts_to_lns = "pts_to_lns"
    arcpy.CopyFeatures_management(joined_contours, "joined_cont_copy")
    arcpy.MakeFeatureLayer_management(end_pts_from_lns, "end_pts_from_lns")
    for i in fids_list:
        # для выборки необходимы круглые скобки, заменяем
        pts_ids = "OBJECTID IN " + str([i[0], i[1]]).replace('[','(').replace(']',')')
        selected_pairs = arcpy.SelectLayerByAttribute_management("end_pts_from_lns", "NEW_SELECTION", pts_ids)
        arcpy.PointsToLine_management(selected_pairs, "pts_to_lns")
        # создание поля для записи высоты
        arcpy.AddField_management(pts_to_lns, "Habs", "FLOAT")
        addH = arcpy.da.UpdateCursor("pts_to_lns", ["OBJECTID", "Shape", "Shape_Length", "Habs"])
        # запись высоты в поле Habs
        for h in addH:
            h[3] = i[2]
            addH.updateRow(h)
        del h, addH

        # присоединение каждой полилинии к основному массиву
        arcpy.Append_management(pts_to_lns, "joined_cont_copy", "NO_TEST")

    arcpy.AddMessage("Join input supplementary lines with segments of middle lines")
    # объединение исходных горизонталей с фрагментами по значению высоты (слияние по атрибуту)
    dissolved_final_lns = "dissolved_final_lns"
    arcpy.Dissolve_management("joined_cont_copy", dissolved_final_lns, "Habs")

    # перевод изолиний в отдельные объекты
    final_lines = output_contours
    arcpy.MultipartToSinglepart_management(dissolved_final_lns, final_lines)

    # удаление ненужных колонок с атрибутами
    arcpy.DeleteField_management(output_contours, "ORIG_FID")

    # создание shape-файла с новыми горизонталями
    # arcpy.FeatureClassToFeatureClass_conversion(final_lines, , "finally_processed_lns.shp") 


if __name__ == '__main__':
    # разделение выполнения операции на несколько процессов для ускорения работы
    arcpy.env.parallelProcessingFactor = "0"
    # активация возможности перезаписи файлов
    arcpy.env.overwriteOutput = True
    args = tuple(arcpy.GetParameterAsText(i) for i in range(arcpy.GetArgumentCount()))
    contour_prolongation(*args)
