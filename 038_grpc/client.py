import converter_pb2
import converter_pb2_grpc
import grpc
import streamlit as st


def send_message(message: str):
    with grpc.insecure_channel("grpc-server:50051") as channel:
        stub = converter_pb2_grpc.ConverterStub(channel)
        request = converter_pb2.ConvertRequest(message=message)
        response = stub.Convert(request)
        return response.result


st.title("Simple gRPC Demo")

message = st.text_input("Enter message:")
if st.button("Convert to Uppercase"):
    result = send_message(message)
    st.write("Result:", result)
