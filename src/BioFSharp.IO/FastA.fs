﻿namespace BioFSharp.IO

open System
open FSharpAux
open FSharpAux.IO

module FastA =
    open System.IO
            
    /// Fasta item contains header and sequence
    type FastaItem<'a> = {
        Header    : string;
        Sequence  : 'a;       
    }

    let toTaggedSequence (fsa:FastaItem<'S>) =
        BioFSharp.TaggedSequence.create fsa.Header fsa.Sequence
            
    /// Creates with header line and sequence.
    let createFastaItem header sequence =
        { Header = header; Sequence = sequence }

        
    // Conditon of grouping lines
    let private same_group l =             
        not (String.length l = 0 || l.[0] <> '>')
    

    
    /// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
    let fromFileEnumerator (converter:seq<char>-> 'a) (fileEnumerator) =
        // Matches grouped lines and concatenates them
        let record d (converter:seq<char>-> 'a) = 
            match d with
            | [] -> raise (System.Exception "Incorrect FASTA format")
            | (h:string) :: t when h.StartsWith ">" ->  let header = h .Remove(0,1)
                                                        let sequence = (Seq.concat t) |> converter
                                                        createFastaItem header sequence
                                                        
            | h :: _ -> raise (System.Exception "Incorrect FASTA format")        
    
        // main
        fileEnumerator
        |> Seq.filter (fun (l:string) -> not (l.StartsWith ";" || l.StartsWith "#"))
        |> Seq.groupWhen same_group 
        |> Seq.map (fun l -> record (List.ofSeq l) converter)
    

    /// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
    let fromFile converter (filePath) =
        FileIO.readFile filePath
        |> fromFileEnumerator converter



    /// Reads FastaItem from gzFile. Converter determines type of sequence by converting seq<char> -> type
    let fromGzipFile converter (filePath) =
        FileIO.readFileGZip filePath
        |> fromFileEnumerator converter


    /// Writes FastaItem to stream. Converter determines type of sequence by converting type -> char
    /// The passed stream stays open and is not disposed after writing to it.
    /// If you want to reuse the stream (e.g. you are not writing to a file stream but a memory stream that gets used afterwards)
    /// you have to reset the position with `stream.Seek(0L, SeekOrigin.Begin)`
    let writeToStream (toString:'T -> char) (stream:Stream) (data:seq<FastaItem<#seq<'T>>>) =
        let toChunks (w:System.IO.StreamWriter) (length:int) (source: seq<'T>) =    
            use ie = source.GetEnumerator()
            let sourceIsEmpty = ref false
            let builder = System.Text.StringBuilder(length)
            let rec loop () =        
                    if ie.MoveNext () then                
                        builder.Append(toString ie.Current) |> ignore
                        for x in 2 .. length do
                            if ie.MoveNext() then
                                builder.Append(toString ie.Current) |> ignore
                            else
                                sourceIsEmpty := true                
                
                        match !sourceIsEmpty with
                        | false -> // writer builder
                                   w.WriteLine(builder.ToString())
                                   builder.Clear() |> ignore
                                   loop ()
                        | true  -> w.WriteLine(builder.ToString())
                                   ()
                    else
                        w.Flush()
        
            loop ()
        use sWriter = new System.IO.StreamWriter(stream,Text.UTF8Encoding(false,true),4096,true)
        data
        |> Seq.iter (fun (i:FastaItem<_>) ->
            sWriter.WriteLine(">" + i.Header)
            toChunks sWriter 80 i.Sequence
        ) 


    /// Writes FastaItem to file. Converter determines type of sequence by converting type -> char. If file already exists the data is overwritten.
    let write (toString:'T -> char) (filePath:string) (data:seq<FastaItem<#seq<'T>>>) =
        let file = new FileStream(filePath,FileMode.Create)
        writeToStream toString file data
        file.Dispose()

    /// Writes FastaItem to file. Converter determines type of sequence by converting type -> char. If file already exists the data is appended.
    let writeAndAppend (toString:'T -> char) (filePath:string) (data:seq<FastaItem<#seq<'T>>>) =
        let file = new FileStream(filePath,FileMode.Append)
        writeToStream toString file data   
        file.Dispose()

    /// Converts FastaItem to string. Converter determines type of sequence by converting type -> char
    let toString (toString:'T -> char) (data:seq<FastaItem<#seq<'T>>>) =
        let toChunks (length:int) (source: seq<'T>) (head:string)=    
            let ie = source.GetEnumerator()
            let sourceIsEmpty = ref false
            let builder = System.Text.StringBuilder(length)        
            seq {
                yield sprintf ">%s" head
                while ie.MoveNext () do                
                            builder.Append(toString ie.Current) |> ignore
                            for x in 2 .. length do
                                if ie.MoveNext() then
                                    builder.Append(toString ie.Current) |> ignore
                                else
                                    sourceIsEmpty := true                
            
                            match !sourceIsEmpty with
                            | false -> // writer builder
                                        yield (builder.ToString())
                                        builder.Clear() |> ignore
                            | true  -> yield (builder.ToString())
                }
        data
        |> Seq.map (fun (i:FastaItem<_>) -> toChunks 80 i.Sequence i.Header)
        |> Seq.concat

